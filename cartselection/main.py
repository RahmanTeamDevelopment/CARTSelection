
import helper
import sys


def main(options):
    """Main function"""

    # Initialize output files
    selected, missing, log = helper.initialize_output_files(options)

    # Read input datasets
    gene_dict, appris, appris_principal, refseq_db, excluded, refseqscan = helper.readInputData(options)

    # Process genes
    process_genes(options, gene_dict, appris, appris_principal, refseq_db, refseqscan, selected, missing, log)


def process_genes(options, gene_dict, appris, appris_principal, refseq_db, refseqscan, selected, missing, log):
    """Iterate through inputed genes and process each one"""

    sys.stdout.write('Processing genes ... ')
    sys.stdout.flush()

    cartidx = 10000
    # Iterate through HGNC IDs from input file
    for hgncid in open(options.hgnc):

        hgncid = hgncid.strip()
        if hgncid == '': continue

        # Log gene name
        log.write('-' * 90 + '\n')
        log.write('\n\n' + '=' * 90 + '\n')
        log.write('Gene: ' + hgncid + '\n')

        # Get APPRIS PRINCIPAL isoforms
        nms = get_appris_principal_isoforms(gene_dict, appris, hgncid, missing, log)
        if nms is None:
            continue

        # Get transcript data from RefSeq db
        transcripts = get_transcript_data_from_refseq_db(refseq_db, nms, missing, log, hgncid)
        if transcripts is None:
            continue

        # CART series and genome build
        cartseries = options.series.upper()
        build = refseq_db.build[-2:] if options.build is None else options.build

        # Select and output transcript
        cartidx = select_transcript(transcripts, refseqscan, appris_principal, selected, missing, log, hgncid, build, cartseries, cartidx)

    # Close output files
    selected.close()
    missing.close()
    log.close()


def get_appris_principal_isoforms(gene_dict, appris, hgncid, missing, log):
    """Get APPRIS PRINCIPAL isoforms"""

    # Translate HGNC ID to NCBI gene ID
    ncbi_geneid = helper.translate_gene_id(gene_dict, hgncid, missing, log)
    if ncbi_geneid is None:
        return

    # Get APPRIS principal isoforms for the gene
    appris_principal = helper.get_appris_principal_isoforms_for_gene(appris, ncbi_geneid, missing, log, hgncid)
    if appris_principal is None:
        return

    # Check if only XMs are present in APPRIS
    nms = helper.get_nms(appris_principal, missing, log, hgncid)
    if nms is None:
        return
    return nms


def get_transcript_data_from_refseq_db(refseq_db, nms, missing, log, hgncid):
    """Get transcript data from RefSeq db"""

    # Get NMs that are included in the RefSeq DB
    nms_in_refseq = helper.check_nms_in_refseq_db(refseq_db, nms, missing, log, hgncid)
    if nms_in_refseq is None:
        return

    # Create Transcript objects
    transcripts = helper.create_transcript_objects(refseq_db, nms_in_refseq)

    # Get transcripts with consistent HGNC ID
    transcripts = helper.get_transcripts_with_consistent_hgncid(transcripts, hgncid, missing, log)
    if transcripts is None:
        return
    return transcripts


def select_transcript(transcripts, refseqscan, appris_principal, selected, missing, log, hgncid, build, cartseries, cartidx):
    """Select transcript"""

    # Single transcript to select
    if len(transcripts) == 1:
        cartidx += 1
        helper.output_single_transcript(transcripts, refseqscan, appris_principal, selected, log, hgncid, build, cartseries, cartidx)
        return cartidx

    # Multiple transcripts to select from
    elif len(transcripts) > 1:
        return helper.select_from_multiple_candidates(transcripts, refseqscan, appris_principal, selected, missing, log, hgncid, build, cartseries, cartidx)

