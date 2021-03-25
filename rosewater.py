from pathlib import Path
import pandas as pd
import pybedtools as pyb
from pybedtools.featurefuncs import TSS
import click


def parse_annotation(anno_gtf, out_dir):

    # Get base filename.
    base = Path(anno_gtf).stem

    anno = pyb.BedTool(anno_gtf)
    all_tss = anno.filter(lambda b: b[2] == "transcript").each(
        TSS, upstream=0, downstream=1).saveas()
    pc_tss = all_tss.filter(
        lambda x: x.attrs["transcript_type"] == "protein_coding").saveas()

    return all_tss, pc_tss


def filter_annotations(sample, tpm_threshold, tpm_file, all_tss, pc_tss):

    # Get and filter TPMs for sample.
    tpms = pd.read_table(tpm_file)
    tpms = dict(zip(tpms["gene"], tpms[sample]))
    filt_tpms = [key for (key, value) in tpms.items()
                 if value >= tpm_threshold]

    filt_all_tss = all_tss.filter(
        lambda x: x.attrs['gene_name'] in filt_tpms).sort().saveas()
    filt_pc_tss = pc_tss.filter(
        lambda x: x.attrs['gene_name'] in filt_tpms).sort().saveas()

    return tpms, filt_all_tss, filt_pc_tss


def gene_intersects(enh_file, filt_all_tss, filt_pc_tss):
    enh_df = pd.read_table(enh_file)
    enh_df.drop(enh_df.columns[0], axis=1, inplace=True)
    enhancers = enh_df.to_string(header=False, index=False, columns=[
                                 "CHROM", "START", "STOP", "CLOSEST_GENE", "enhancerRank"])
    enhancers = pyb.BedTool(enhancers, from_string=True).sort().saveas()

    ovlps = enhancers.intersect(filt_all_tss, loj=True).each(parse_ovlp_gene_info).merge(c=[4, 5, 15, 16, 17],
                                                                                         o=["distinct", "distinct", "distinct", "collapse", "distinct"]).saveas()

    return ovlps


def parse_ovlp_gene_info(feature):
    # If TSS overlaps enhancer, get relevant transcript and gene info.
    fields = feature.fields
    info = fields[-1]
    symb = "."
    t_id = "."
    ttype = "."
    if info != ".":
        symb = info.strip().split(";")[4]
        symb = symb.split()[1]
        symb = symb.strip('"')

        t_id = info.strip().split(";")[1]
        t_id = t_id.split()[1]
        t_id = t_id.strip('"')
        ttype = info.strip().split(";")[5]
        ttype = ttype.split()[1]
        ttype = ttype.strip('"')

    fields.extend([t_id, ttype, symb])
    interval = pyb.cbedtools.create_interval_from_list(fields)
    return interval


def get_closest_genes(ovlps, filt_all_tss, filt_pc_tss):
    closest_genes = ovlps.closest(filt_all_tss, d=True, t="first").closest(filt_pc_tss, d=True, t="first").each(
        parse_closest_gene_info).cut([0, 1, 2, 4, 3, 5, 6, 7, 28, 29, 17, 30, 27]).saveas()
    return closest_genes


def parse_closest_gene_info(feature):
    # For nearest all and protein coding TSS, get gene type and symbol.
    fields = feature.fields
    all_info = fields[16]
    pc_info = fields[-2]

    all_ttype = all_info.strip().split(";")[5]
    all_ttype = all_ttype.split()[1]
    all_ttype = all_ttype.strip('"')
    all_symb = all_info.strip().split(";")[4]
    all_symb = all_symb.split()[1]
    all_symb = all_symb.strip('"')
    pc_symb = pc_info.strip().split(";")[4]
    pc_symb = pc_symb.split()[1]
    pc_symb = pc_symb.strip('"')

    fields.extend([all_symb, all_ttype, pc_symb])
    interval = pyb.cbedtools.create_interval_from_list(fields)
    return interval


def assign_genes(enh_bedtools, tpms, transcript_dict, genetype_dict, out_file, log_file):

    # For logging.
    total_enh = 0
    no_tss = 0
    one_tss = 0
    two_tss = 0
    three_or_more_tss = 0
    all_diff_from_rose = 0
    pc_diff_from_rose = 0

    # Parse enhancer and get TPM, transcript_type, and gene_type values.
    for enh in enh_bedtools:
        total_enh += 1

        pos = "\t".join(enh.fields[0:4])
        rose_assignment = transcript_dict[enh.fields[4]]
        rose_gtype = genetype_dict[rose_assignment]
        if rose_assignment not in tpms.keys():
            rose_tpm = "NA"
        else:
            rose_tpm = tpms[rose_assignment]
        ovlp_tss = enh.fields[5]
        ovlp_ttypes = enh.fields[6]
        ovlp_genes = enh.fields[7].split(",")
        ovlp_gtypes = "."
        ovlp_tpms = "."
        closest_tss = enh.fields[8]
        closest_tss_ttype = enh.fields[9]
        closest_tss_gtype = genetype_dict[closest_tss]
        if closest_tss not in tpms.keys():
            closest_tss_tpm = "NA"
        else:
            closest_tss_tpm = tpms[closest_tss]

        closest_tss_distance = enh.fields[10]
        closest_pc_tss = enh.fields[11]

        if closest_pc_tss not in tpms.keys():
            closest_pc_tss_tpm = "NA"
        else:
            closest_pc_tss_tpm = tpms[closest_pc_tss]

        closest_pc_tss_distance = enh.fields[-1]

        if ovlp_genes[0] != ".":
            ovlp_gtypes = [genetype_dict[x] for x in ovlp_genes]
            ovlp_tpms = [tpms[x] for x in ovlp_genes]

        # Check for no overlapping TSSs.
        if ovlp_genes[0] == ".":
            no_tss += 1

            final_pc_assignment = closest_pc_tss
            final_pc_assignment_tpm = closest_pc_tss_tpm
            final_all_assignment = closest_tss
            final_all_assignment_gtype = closest_tss_gtype
            final_all_assignment_tpm = closest_tss_tpm

        # If only one TSS is overlapped and it's protein coding, assignment will be the same.
        elif len(ovlp_genes) == 1 and ovlp_gtypes[0] != ".":
            one_tss += 1
            final_all_assignment = ovlp_genes[0]
            final_all_assignment_gtype = ovlp_gtypes[0]
            final_all_assignment_tpm = ovlp_tpms[0]

            if ovlp_gtypes[0] == "protein_coding":
                final_pc_assignment = ovlp_genes[0]
                final_pc_assignment_tpm = ovlp_tpms[0]
            else:
                final_pc_assignment = closest_pc_tss
                final_pc_assignment_tpm = closest_pc_tss_tpm

        elif len(ovlp_genes) > 1:
            if len(ovlp_genes) == 2:
                two_tss += 1
            elif len(ovlp_genes) > 2:
                three_or_more_tss += 1

            # Get tpms, find max, get associated gene and change assignment if necessary.
            index_max = max(range(len(ovlp_tpms)), key=ovlp_tpms.__getitem__)
            highest_tpm_ovlp = ovlp_genes[index_max]
            highest_tpm_gtype = ovlp_gtypes[index_max]

            if highest_tpm_gtype == "protein_coding":
                final_pc_assignment = highest_tpm_ovlp
                final_pc_assignment_tpm = max(ovlp_tpms)

                final_all_assignment = highest_tpm_ovlp
                final_all_assignment_gtype = highest_tpm_gtype
                final_all_assignment_tpm = max(ovlp_tpms)

            elif highest_tpm_gtype != "protein_coding":

                # Get max protein coding gene.
                pc_gtypes_index = [idx for idx, element in enumerate(
                    ovlp_gtypes) if element == "protein_coding"]
                pc_ovlp_tpms = [ovlp_tpms[x] for x in pc_gtypes_index]
                pc_gene = [ovlp_genes[x] for x in pc_gtypes_index]
                pc_index_max = max(range(len(pc_ovlp_tpms)),
                                   key=pc_ovlp_tpms.__getitem__)
                highest_pc_tpm = pc_ovlp_tpms[pc_index_max]
                highest_pc_tpm_ovlp = pc_gene[pc_index_max]

                final_all_assignment = highest_tpm_ovlp
                final_all_assignment_gtype = highest_tpm_gtype
                final_all_assignment_tpm = max(ovlp_tpms)

                if highest_pc_tpm > closest_pc_tss_tpm:
                    final_pc_assignment = highest_pc_tpm_ovlp
                    final_pc_assignment_tpm = highest_pc_tpm
                else:
                    final_pc_assignment = closest_pc_tss
                    final_pc_assignment_tpm = closest_pc_tss_tpm

        # Construct and print new line.
        print("\t".join([pos, rose_assignment, rose_gtype, str(rose_tpm),
                         ovlp_tss, ovlp_ttypes, ",".join(ovlp_genes), ",".join(
                             ovlp_gtypes), ",".join([str(x) for x in ovlp_tpms]),
                         closest_tss, closest_tss_ttype, closest_tss_gtype, str(
                             closest_tss_tpm), str(closest_tss_distance),
                         closest_pc_tss, str(closest_pc_tss_tpm), str(
                             closest_pc_tss_distance),
                         final_all_assignment, final_all_assignment_gtype, str(final_all_assignment_tpm), 
                         final_pc_assignment, str(final_pc_assignment_tpm)]),
              file=out_file)

        if final_all_assignment != rose_assignment:
            all_diff_from_rose += 1
        if final_pc_assignment != rose_assignment:
            pc_diff_from_rose += 1

        # Print stats to log.
    print("Total enhancers: " + str(total_enh), file=log_file)
    print("Average enhancer size (bp): " + str(enh_bedtools.total_coverage() /
                                         enh_bedtools.count()), file=log_file)
    print("Enhancers overlapping no TSS: " + str(no_tss) + ", " +
          str(100 * (no_tss / total_enh)) + "% of total enhancers", file=log_file)
    print("Enhancers overlapping one TSS: " + str(one_tss) + ", " +
          str(100 * (one_tss / total_enh)) + "% of total enhancers", file=log_file)
    print("Enhancers overlapping two TSSes: " + str(two_tss) + ", " +
          str(100 * (two_tss / total_enh)) + "% of total enhancers", file=log_file)
    print("Enhancers overlapping three or more TSSes: " + str(three_or_more_tss) + ", " + str(100 * (three_or_more_tss / total_enh)) + "% of total enhancers",
          file=log_file)
    print("Enhancers using all TSSes where assignment differed from ROSE: " + str(all_diff_from_rose) + ", " +
          str(100 * (all_diff_from_rose / total_enh)) + "% of total enhancers", file=log_file)
    print("Enhancers using only protein-coding TSSes where assignment differed from ROSE: " + str(pc_diff_from_rose) + ", " +
          str(100 * (pc_diff_from_rose / total_enh)) + "% of total enhancers", file=log_file)


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-a", "--annotation", help="Path to annotation GTF file.", required=True, type=click.Path())
@click.option("-s", "--sample", help="Sample name that should match a column in the TPM file.", required=True, type=str)
@click.option("-e", "--enh_file", help="Output from ROSE ending with 'ENHANCER_TO_GENE.txt'.", required=True, type=click.Path())
@click.option("-t", "--tpm_file", help="A file containing a matrix of TPMs with genes as rows and samples as columns. The gene label column should be named 'gene'.",
              required=True, type=click.Path())
@click.option("-th", "--tpm_threshold", default=5, help="The minimum TPM to retain genes for assignment.", show_default=True, type=float)
@click.option("-o", "--output_dir", default="./EnhancerGeneAssignments", help="The output directory.", show_default=True, type=click.Path())
@click.version_option()
def rosewater(annotation, sample, enh_file, tpm_file, tpm_threshold, output_dir):
    """rosewater assigns genes to ROSE output in an expression-aware manner."""

    # Make output directory and open file for logging.
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    base = Path(enh_file).stem
    log_file = open(
        Path(output_dir, base + ".rosewater.GeneAssignment.log"), "w")

    print("Annotation: " + annotation, file=log_file)
    print("Sample: " + sample, file=log_file)
    print("ROSE file: " + enh_file, file=log_file)
    print("TPM file: " + tpm_file, file=log_file)
    print("TPM threshold: " + str(tpm_threshold), file=log_file)
    print("Output directory: " + output_dir + "\n", file=log_file)

    # Get all TSSes and only protein coding ones.
    all_tss, pc_tss = parse_annotation(annotation, output_dir)
    print("Total TSSes: " + str(len(all_tss)), file=log_file)
    print("Protein coding TSSes: " + str(len(pc_tss)), file=log_file)

    # Filter them to get only genes that meet the TPM threshold.
    tpms, filt_all_tss, filt_pc_tss = filter_annotations(
        sample, tpm_threshold, tpm_file, all_tss, pc_tss)
    print("TSSes that met TPM threshold (" +
          str(tpm_threshold) + "): " + str(len(filt_all_tss)), file=log_file)
    print("Protein coding TSSes that met TPM threshold (" +
          str(tpm_threshold) + "): " + str(len(filt_pc_tss)) + "\n", file=log_file)

    # Make dicts of transcript IDs and gene symbols and gene types.
    transcript_dict = {transcript.attrs["transcript_id"]: transcript.attrs["gene_name"] for transcript in all_tss}
    genetype_dict = {transcript.attrs["gene_name"]: transcript.attrs["gene_type"] for transcript in all_tss}

    # Get overlapping TSSes.
    ovlps = gene_intersects(enh_file, filt_all_tss, filt_pc_tss)

    # Get closest genes with all TSSes and only protein-coding ones.
    closest_genes = get_closest_genes(ovlps, filt_all_tss, filt_pc_tss)

    # Collect and print output for each SE.
    header = ["chrom", "start", "stop", "enhancerRank", "ROSE_assignment", "ROSE_assignment_genetype", "ROSE_assignment_TPM",
              "overlapping_TSS", "overlapping_TSS_transcripttype", "overlapping_genes", "overlapping_genetypes", "overlapping_genes_TPMs",
              "closest_TSS", "closest_TSS_transcripttype", "closest_TSS_genetype", "closest_TSS_gene_TPM", "closest_TSS_distance",
              "closest_proteincoding_TSS", "closest_proteincoding_TSS_gene_TPM", "closest_proteincoding_TSS_distance",
              "final_allgene_assignment", "final_allgene_assignment_genetype","final_allgene_assignment_TPM", "final_proteincoding_assignment", "final_proteincoding_assignment_TPM"]
    out_file = open(
        Path(output_dir, base + ".rosewater.GeneAssignment.bed"), "w")
    print("\t".join(header), file=out_file)

    assign_genes(closest_genes, tpms, transcript_dict,
                 genetype_dict, out_file, log_file)

    out_file.close()
    log_file.close()
