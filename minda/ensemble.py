import sys
from ast import parse
from collections import Counter
from datetime import datetime
import pandas as pd
import numpy as np
import re
import gzip
from minda.decompose import _is_vcf_gz

def parse_bnd_alt(alt_string):
    '''
    Parse the BND alt string and return separators and region
    adapted from svtools package
    '''
    # NOTE The below is ugly but intended to match things like [2:222[ and capture the brackets
    result = re.findall(r'([][])(.+?)([][])', alt_string)
    assert result, "%s\n" % alt_string
    sep1, _ , sep2 = result[0]
    assert sep1 == sep2
    return sep1

def _infer_strands(svtype, alt):
    """
    infer SV strands from ALT
    """
    valid_svtypes = {"DEL", "INS", "DUP", "INV", "BND"}
    assert svtype in valid_svtypes, f"invalid svtype: {svtype}. must be one of {valid_svtypes}"

    orientation1 = orientation2 = "+"
    if svtype in ("DEL", "INS") or alt in ("<DEL>", "<INS>"):
        orientation2 = "-"
    elif svtype == "DUP" or "DUP" in alt:
        orientation1 = "-"
    else:
        sep = parse_bnd_alt(alt)
        if alt.startswith(sep):
            orientation1 = "-"
        if sep == "[":
            orientation2 = "-"
    return orientation1+orientation2

def _add_columns(ensemble_df, vaf):
    # create a column of list of prefixed IDs for each locus group
    key_columns = ['locus_group_x','locus_group_y']
    value_columns = ['ID_x', 'ID_y']
    column_suffixes = ['x','y']
    
    for i in range(len(key_columns)):
        locus_group = key_columns[i]
        id = value_columns[i]
        column_suffix = column_suffixes[i]
        
        keys =  ensemble_df[f'{locus_group}'].to_list()
        values = ensemble_df[f'{id}'].to_list()
        minda_values = ensemble_df['Minda_ID'].to_list()
        caller_names =  ensemble_df.caller_names.to_list()
        
        id_dict = {}
        for key, value in zip(keys, values):
            if key not in id_dict:
                id_dict[key] = []
            id_dict[key].append(value)
    
        minda_id_dict = {}
        for key, value in zip(keys, minda_values):
            if key not in minda_id_dict:
                minda_id_dict[key] = []
            minda_id_dict[key].append(value)

        ensemble_df[f'ID_list_{column_suffix}'] = ensemble_df[locus_group].map(id_dict)
        ensemble_df[f'Minda_ID_list_{column_suffix}'] = ensemble_df[locus_group].map(minda_id_dict)

        
    # create dict for SV type
    values = ensemble_df.SVTYPE.to_list()
    svtype_dict = {}
    for key, value in zip(keys, values):
        if key not in svtype_dict:
            svtype_dict[key] = []
        svtype_dict[key].append(value)
        
    most_common_svtpye_dict = {k:Counter(v).most_common(1)[0][0] for (k,v) in svtype_dict.items()}
    ensemble_df['SVTYPE'] = ensemble_df['locus_group_y'].map(most_common_svtpye_dict)

    # create dict of  vafs
    if vaf != None:
        values = ensemble_df.VAF.to_list()
        vaf_dict = {}
        for key, value in zip(keys, values):
            if key not in vaf_dict:
                vaf_dict[key] = []
            vaf_dict[key].append(value)
    
        ensemble_df['VAFs'] = ensemble_df['locus_group_y'].map(vaf_dict)
    
    return ensemble_df


def _get_ensemble_df(decomposed_dfs_list, caller_names, tolerance, vaf, out_dir, sample_name, args, multimatch):
    
    dfs_1 = [dfs_list[0] for dfs_list in decomposed_dfs_list]
    dfs_2 = [dfs_list[1] for dfs_list in decomposed_dfs_list]
    dfs_list = [dfs_1, dfs_2]
    
    # create stat dfs
    start_dfs_list = []
    start_dfs = pd.concat(dfs_1).reset_index(drop=True)
    start_dfs = start_dfs[['#CHROM', 'POS', 'ID', 'Minda_ID', 'INFO', 'SVTYPE', 'SVLEN', 'REF', 'ALT']].sort_values(['#CHROM', 'POS'])
    
    start_dfs['diff_x'] = start_dfs.groupby('#CHROM').POS.diff().fillna(9999)
    diffs = start_dfs['diff_x'].to_list()

    # group start loci
    loci = []
    count = 1
    for diff in diffs:
        if diff >= tolerance:
            locus = count
            loci.append(locus)
            count += 1
        else:
            locus = count - 1 
            loci.append(locus)
    
    
    start_dfs['locus_group_x'] = loci
    start_dfs['median'] = start_dfs.groupby('locus_group_x')['POS'].transform('median').astype('int')

    # create end dfs
    end_dfs = pd.concat(dfs_2).reset_index(drop=True)

    #ensemble_df = start_dfs.merge(end_dfs, on=['SVTYPE', 'SVLEN','Minda_ID'])
    ensemble_df = start_dfs.merge(end_dfs, on=['SVTYPE', 'Minda_ID'])
    ensemble_df = ensemble_df.sort_values(['locus_group_x','#CHROM_y', 'POS_y'])
    ensemble_df ['diff_y'] = ensemble_df.groupby(['locus_group_x','#CHROM_y']).POS_y.diff().abs().fillna(9999)
    diffs = ensemble_df['diff_y'].to_list()
    caller_names = ensemble_df['Minda_ID'].apply(lambda x: x.rsplit('_', 1)[0]).tolist()
    ensemble_df['caller_names']= caller_names

    # group end loci
    locus_callers = []
    loci = []
    count = 1
    for i in range(len(diffs)):
        diff = diffs[i]
        caller_name = caller_names[i]
        
        #create new end locus (sublocus 1)
        if diff >= tolerance: 
            locus_callers.clear()
            locus_callers.append(caller_name)
            locus = str(count) + "_1"
            loci.append(locus)
            count += 1
            
        # add new sublocus (sublocus determined by how many calls the caller makes at a given locus)
        elif multimatch == False and diff < tolerance and caller_name in locus_callers:
            #elif diff < tolerance and caller_name in locus_callers: 
                locus_callers.append(caller_name)
                sub_group = locus_callers.count(caller_name)
                locus = str(count-1) + "_" + str(sub_group)
                loci.append(locus)
            
        # add to existing sublocus 1   
        else: 
            locus = str(count - 1) + "_1" 
            loci.append(locus)
            locus_callers.append(caller_name)
    
    # add locus_group_y, median POS, caller ID, Minda ID, SV type, VAF columns
    ensemble_df['locus_group_y'] = loci
    ensemble_df['median'] = ensemble_df.groupby('locus_group_y')['POS_y'].transform('median').astype('int')
    ensemble_df = _add_columns(ensemble_df, vaf)
    if vaf != None:
        ensemble_df['VAF'] = ensemble_df.groupby('locus_group_y')['VAF'].transform('median')
    else:
        ensemble_df['VAF'] = np.nan
   
    ensemble_df = ensemble_df.drop_duplicates(['locus_group_x', 'locus_group_y']).reset_index(drop=True)
    
    return ensemble_df


def _get_ensemble_call_column(support_df, conditions):
    column_names = []
    condition_count = 0
    condition_columns = []
    query_list = []
    for i in range(len(conditions)):
        
        if i % 2 == 0:
            operator = conditions[i][1]
            number = str(conditions[i][2])
            
            nested_caller_columns = conditions[i][0]
            nested_type = type(nested_caller_columns[0])
            
            condition = chr(ord('A') + condition_count)
            column_name = f'condition_{condition}'
    
            nested_columns_count = 1
            sub_condition_columns = []
            if nested_type == list:
                for j in range(len(nested_caller_columns)):
                    caller_columns = nested_caller_columns[j]
                    sub_column_name = f'condition_{nested_columns_count}_{condition}'
                    support_df[f'{sub_column_name}'] = support_df[caller_columns].any(axis=1)
                    nested_columns_count += 1
                    sub_condition_columns.append(sub_column_name)
                support_df[f'{column_name}'] = support_df[sub_condition_columns].sum(axis=1)
                
            else:
                support_df[f'{column_name}'] = support_df[nested_caller_columns].sum(axis=1)
            condition_count += 1
            condition_columns.append(column_name)
            query_list.extend([column_name, operator, number])
        else:
            query_list.extend(conditions[i])
    
    query = ' '.join(query_list)
    mask = support_df.eval(query) 
    #support_df['ensemble'] = mask 
    support_df.insert(loc=12, column='ensemble', value=mask)
    return support_df

def _get_contigs(vcf_list):
    contig_dict = {}
    for vcf in vcf_list:
        is_vcf_gz = _is_vcf_gz(vcf)
        open_func = gzip.open if is_vcf_gz else open
        mode = 'rt' if is_vcf_gz else 'r'

        with open_func(vcf, mode) as file:
            for line in file:
                if not line.startswith("##"):
                    break
                if line.startswith("##contig"):
                    pattern = r'ID=([^,>]+),length=([^,>]+)'
                    id_length_tuple = re.findall(pattern, line)[0]
                    chr_id = id_length_tuple[0]
                    length = int(id_length_tuple[1])
                    if chr_id not in contig_dict:
                        contig_dict[chr_id] = length
                    else:
                        value = contig_dict[chr_id]
                        max_length = max(value, length)
                        if value != max_length:
                            print(f'Contig ID {chr_id} has lengths {value} and {max_length}; {max_length} will be used in ensemble.vcf header.')
                        contig_dict[chr_id] = max_length
    return contig_dict

def _get_ensemble_vcf(vcf_list, support_df, out_dir, sample_name, args, vaf, version):
    vcf_df = support_df[support_df['ensemble'] == True].reset_index(drop=True).copy()
    vcf_df['ID'] = f'Minda_' + (vcf_df.index + 1).astype(str)
    vcf_df['QUAL'] = "."
    vcf_df['FILTER'] = "PASS"

    if vaf != None:
        vcf_df['INFO'] = ['SVLEN=' + str(svlen) + ';SVTYPE=' + svtype + \
                          ';CHR2=' + str(chr2) + ';END=' + str(end) + \
                          ';STRANDS=' + strands + \
                          ';SUPP_VEC=' + ','.join(map(str, supp_vec)) + ';VAF=' + str(vaf) \
                          for svlen, svtype, chr2, end, strands, supp_vec, vaf in zip(vcf_df['SVLEN'],vcf_df['SVTYPE'], vcf_df['#CHROM_y'], vcf_df['POS_y'], vcf_df['STRANDS'], vcf_df['ID_list_y'], vcf_df['VAF'])]
        vcf_df = (vcf_df[['#CHROM_x', 'POS_x', 'ID', 'REF_x', 'ALT_x', 'QUAL', 'FILTER','INFO']]
                  .rename(columns={'#CHROM_x':"#CHROM", "POS_x":"POS", "REF_x":"REF", "ALT_x": "ALT"}))
    else:
        vcf_df['INFO'] = ['SVLEN=' + str(svlen) + ';SVTYPE=' + svtype + \
                          ';CHR2=' + str(chr2) + ';END=' + str(end) + \
                          ';STRANDS=' + strands + \
                          ';SUPP_VEC=' + ','.join(map(str, supp_vec)) \
                          for svlen, svtype, chr2, end, strands, supp_vec in zip(vcf_df['SVLEN'],vcf_df['SVTYPE'], vcf_df['#CHROM_y'], vcf_df['POS_y'], vcf_df['STRANDS'], vcf_df['ID_list_y'])]
        vcf_df = (vcf_df[['#CHROM_x', 'POS_x', 'ID', 'REF_x', 'ALT_x', 'QUAL', 'FILTER','INFO']]
                  .rename(columns={'#CHROM_x':"#CHROM", "POS_x":"POS", "REF_x": "REF", "ALT_x": "ALT"}))
    date = datetime.today().strftime('%Y-%m-%d')
    with open(f'{out_dir}/{sample_name}_minda_ensemble.vcf', 'w') as file:
        file.write(f'##fileformat=VCFv4.2\n##fileDate={date}\n##source=MindaV{version}\n')
        command_str = " ".join(sys.argv)
        file.write(f"##CommandLine= {command_str}\n")
        contig_dict = _get_contigs(vcf_list)
        for key, value in contig_dict.items():
            file.write(f'##contig=<ID={key},length={value}>\n')
        file.write('##ALT=<ID=DEL,Description="Deletion">\n##ALT=<ID=INS,Description="Insertion">\n##ALT=<ID=DUP,Description="Duplication">\n##ALT=<ID=INV,Description="Inversion">\n')
        file.write('##FILTER=<ID=PASS,Description="Default">\n')
        file.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
                   '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">\n'
                   '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">\n'
                   '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
                   '##INFO=<ID=STRANDS,Number=1,Type=String,Description="Breakpoint strandedness (first_to_second)">\n'
                   '##INFO=<ID=SUPP_VEC,Number=.,Type=String,Description="IDs of support records">\n')
        if vaf != None:
            file.write('##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency">\n')
        vcf_df.to_csv(file, sep="\t", index=False)


def get_support_df(vcf_list, decomposed_dfs_list, caller_names, tolerance, conditions, vaf, command, out_dir, sample_name, args, version, multimatch):
    ensemble_df = _get_ensemble_df(decomposed_dfs_list, caller_names, tolerance, vaf, out_dir, sample_name, args, multimatch)
    
    minda_id_x_lists = ensemble_df.Minda_ID_list_x.to_list()
    minda_id_y_lists = ensemble_df.Minda_ID_list_y.to_list()
    
    # check that both start & end have same IDs
    for caller_name in caller_names:
        caller_column = []
        for i in range(len(minda_id_x_lists)):
            minda_id_x_list = minda_id_x_lists[i]
            minda_id_y_list = minda_id_y_lists[i]
            intersect_list = list(set(minda_id_x_list).intersection(set(minda_id_y_list)))
            call_boolean  = any(value.startswith(caller_name) for value in intersect_list)
            caller_column.append(call_boolean)
        ensemble_df[f'{caller_name}'] = caller_column
    # add in a column for strands
    strands = []
    for row in ensemble_df.itertuples():
        alt = row.ALT_x
        svtype = row.SVTYPE
        strands1 = _infer_strands(svtype, alt)
        strands.append(strands1)
    ensemble_df["STRANDS"] = strands

    column_names = ['#CHROM_x', 'POS_x', 'locus_group_x', 'ID_list_x',
                    '#CHROM_y', 'POS_y', 'locus_group_y', 'ID_list_y',
                    'REF_x', 'REF_y', 'ALT_x', 'ALT_y', 'STRANDS',
                    'SVTYPE', 'SVLEN', 'VAF', 'Minda_ID_list_y'] + caller_names
    
    support_df = ensemble_df[column_names].rename(columns={"Minda_ID_list_y": "Minda_IDs"}).copy()
    #if command == "ensemble":
    support_df = _get_ensemble_call_column(support_df, conditions)

    # create ensemble vcf
    _get_ensemble_vcf(vcf_list, support_df, out_dir, sample_name, args, vaf, version)

    # create support csv
    support_ex_df = support_df
    support_df.to_csv(f'{out_dir}/{sample_name}_support.tsv', sep='\t', index=False)
    
    return support_df


def add_vaf(row,df,caller_name):
    for item in row['Minda_IDs']:
        if item in df['Minda_ID'].values:
            return df[df['Minda_ID'] == item]['VAF'].values[0]
    return row[f'{caller_name}']
