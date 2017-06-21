import sys


def utr_selection(transcripts, log):
    """UTR selection function"""

    log.write('\nUTR criteria applied to the following transcripts:\n')
    for t in transcripts:
        log.write('- '+t.id+'\n')
        log.write('  ' + str(t.strand) +'; UTR5 exons: ' + str(len(t.utr5_exons)) +'; UTR5 start: ' + str(t.utr5_start()) + '; UTR5 length: ' + str(t.utr5_exonic_content_length()) + '; UTR3 end: ' + str(t.utr3_end()) + '; UTR3 length: ' + str(t.utr3_exonic_content_length()) + '\n')
    log.write('\nUTR criteria selection steps:\n')

    dtypes = []
    dcrits = []
    candidates = transcripts
    while len(candidates) != 1:
        candidates, dtype, dcrit = select_by_type(candidates, log)
        if candidates is None:
            return None, None, None
        dtypes.append(dtype)
        dcrits.append(dcrit)
        log.write('   Filtered to '+str([t.id for t in candidates])+'\n')

    log.write('\n')
    return candidates[0], ','.join(dtypes), ','.join(dcrits)


def select_by_type(transcripts, log):
    """Filter transcripts depending on different type"""

    # Difference types: UTR5_number and UTR5_boundary
    candidates, dtype, dcrit = analyse_difference_type_utr5_number_or_boundary(transcripts, log)
    if candidates is not None:
        return candidates, dtype, dcrit

    # Difference type: UTR_ends
    candidates, dtype, dcrit = analyse_difference_type_UTR_ends(transcripts, log)
    if candidates is not None:
        return candidates, dtype, dcrit

    # Difference type: UTR3
    return analyse_difference_type_UTR3(transcripts, log)


def analyse_difference_type_utr5_number_or_boundary(transcripts, log):
    """Analyse difference types UTR5 number or boundary"""

    # Check if there is difference in the number of UTR5 exons between any transcripts
    diff_utr5_number = False
    for i in range(1, len(transcripts)):
        if len(transcripts[0].utr5_exons) != len(transcripts[i].utr5_exons):
            diff_utr5_number = True
            break

    # Check if there is difference in the UTR5 exon boundaries between any transcripts
    diff_utr5_boundary = False
    if not diff_utr5_number:
        for i in range(1, len(transcripts)):
            if different_utr5_boundary(transcripts[0], transcripts[i]):
                diff_utr5_boundary = True
                break

    # Filtering
    if diff_utr5_number or diff_utr5_boundary:

        # Difference type
        difftype = 'UTR5_number' if diff_utr5_number else 'UTR5_boundary'

        log.write(' * Difference type: ' + difftype + '\n')

        # Choose the transcript with the largest 5' footprint
        candidates = select_largest_5prime_footprint(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR5 start\n')
            return candidates, difftype, 'UTR5_start'

        # Choose the transcript with the longest 5' UTR exonic content
        candidates = select_longest_utr5_exonic_content(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR5 length\n')
            return candidates, difftype, 'UTR5_length'

        # Choose the transcript with the most 5' UTR boundary
        candidates = select_utr5_first_boundary(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR5 boundary\n')
            return candidates, difftype, 'UTR5_boundary'

        sys.exit('Error 1')

    return None, None, None


def analyse_difference_type_UTR_ends(transcripts, log):
    """Analyse difference type UTR ends"""

    # Check if there is difference in the UTR ends between any transcripts
    diff_UTR_ends = False
    ue = transcripts[0].UTR_ends()
    for i in range(1, len(transcripts)):
        if not transcripts[i].UTR_ends() == ue:
            diff_UTR_ends = True
            break

    # Filtering
    if diff_UTR_ends:

        log.write(' * Difference type: UTR_ends\n')

        # Choose the transcript with the largest 5' footprint
        candidates = select_largest_5prime_footprint(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR5 start\n')
            return candidates, 'UTR_ends', 'UTR5_start'

        # Choose the transcript with the largest 3' footprint
        candidates = select_largest_3prime_footprint(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR3 end\n')
            return candidates, 'UTR_ends', 'UTR3_end'

        # Choose the transcript with the longest 3' UTR exonic content
        candidates = select_longest_utr3_exonic_content(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR3 length\n')
            return candidates, 'UTR_ends', 'UTR3_length'

        # Choose the transcript with the most 3' UTR boundary
        candidates = select_utr3_first_boundary(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR3 boundary\n')
            return candidates, 'UTR_ends', 'UTR3_boundary'

        sys.exit('Error 2')

    return None, None, None


def analyse_difference_type_UTR3(transcripts, log):
    """Analyse difference tyoe UTR3"""

    # Check if there is difference in the UTR3 between any transcripts
    diff_utr3 = False
    for i in range(1, len(transcripts)):
        if different_utr3(transcripts[0], transcripts[i]):
            diff_utr3 = True
            break

    # Filtering
    if diff_utr3:

        log.write(' * Difference type: UTR3\n')

        # Choose the transcript with the largest 3' footprint
        candidates = select_largest_3prime_footprint(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR3 end\n')
            return candidates, 'UTR3', 'UTR3_end'

        # Choose the transcript with the longest 3' UTR exonic content
        candidates = select_longest_utr3_exonic_content(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR3 length\n')
            return candidates, 'UTR3', 'UTR3_length'

        # Choose the transcript with the most 3' UTR boundary
        candidates = select_utr3_first_boundary(transcripts)
        if len(candidates) < len(transcripts):
            log.write('   Applied criteria: UTR3 boundary\n')
            return candidates, 'UTR3', 'UTR3_boundary'

        sys.exit('Error 3')

    return None, None, None

def different_utr5_boundary(transcript1, transcript2):
    """Check if two transcripts have different UTR5 exon boundaries"""

    for i in range(len(transcript1.utr5_exons)):
        exon1 = transcript1.utr5_exons[i]
        exon2 = transcript2.utr5_exons[i]
        if exon1[0] != exon2[0] or exon1[1] != exon2[1]: return True
    return False


def different_utr3(transcript1, transcript2):
    """Check if two transcripts have different UTR3 exon boundaries"""

    if len(transcript1.utr3_exons) != len(transcript2.utr3_exons): return True
    for i in range(len(transcript1.utr3_exons)):
        exon1 = transcript1.utr3_exons[i]
        exon2 = transcript2.utr3_exons[i]
        if exon1[0] != exon2[0] or exon1[1] != exon2[1]: return True
    return False


def select_largest_5prime_footprint(transcripts):
    """Select transcript(s) with largest 5' footprint"""

    ret = [transcripts[0]]
    furthest = transcripts[0].utr5_start()
    for i in range(1, len(transcripts)):
        t = transcripts[i]
        utr5footprint = t.utr5_start()
        if utr5footprint == furthest: ret.append(t)
        elif (t.strand == 1 and utr5footprint < furthest) or (t.strand == -1 and utr5footprint > furthest):
            furthest = utr5footprint
            ret = [t]
    return ret


def select_largest_3prime_footprint(transcripts):
    """Select transcript(s) with largest 3' footprint"""
    ret = [transcripts[0]]
    furthest = transcripts[0].utr3_end()
    for i in range(1, len(transcripts)):
        t = transcripts[i]
        utr3footprint = t.utr3_end()
        if utr3footprint == furthest: ret.append(t)
        elif (t.strand == 1 and utr3footprint > furthest) or (t.strand == -1 and utr3footprint < furthest):
            furthest = utr3footprint
            ret = [t]
    return ret


def select_longest_utr5_exonic_content(transcripts):
    """Select transcript(s) with longest UTR5 exonic content"""

    ret = [transcripts[0]]
    longest = transcripts[0].utr5_exonic_content_length()
    for i in range(1, len(transcripts)):
        t = transcripts[i]
        exonic_content_length = t.utr5_exonic_content_length()
        if exonic_content_length == longest: ret.append(t)
        if exonic_content_length > longest:
            longest = exonic_content_length
            ret = [t]
    return ret


def select_longest_utr3_exonic_content(transcripts):
    """Select transcript(s) with longest UTR3 exonic content"""

    ret = [transcripts[0]]
    longest = transcripts[0].utr3_exonic_content_length()
    for i in range(1, len(transcripts)):
        t = transcripts[i]
        exonic_content_length = t.utr3_exonic_content_length()
        if exonic_content_length == longest: ret.append(t)
        if exonic_content_length > longest:
            longest = exonic_content_length
            ret = [t]
    return ret


def select_utr5_first_boundary(transcripts):
    """Select transcript(s) with most 5' UTR boundary"""

    ret = [transcripts[0]]
    if transcripts[0].strand == 1: most5prime_exons = transcripts[0].utr5_exons
    else: most5prime_exons = [[x[1],x[0]] for x in transcripts[0].utr5_exons]
    for i in range(1, len(transcripts)):
        t = transcripts[i]
        if t.strand == 1: exons = t.utr5_exons
        else: exons = [[x[1],x[0]] for x in t.utr5_exons]
        if most5prime_exons == exons: ret.append(t)
        elif (t.strand == 1 and exons < most5prime_exons) or (t.strand == -1 and exons > most5prime_exons):
            ret = [t]
            most5prime_exons = exons
    return ret


def select_utr3_first_boundary(transcripts):
    """Select transcript(s) with most 3' UTR boundary"""

    ret = [transcripts[0]]
    if transcripts[0].strand == 1: most3prime_exons = [[x[1],x[0]] for x in transcripts[0].utr3_exons[::-1]]
    else: most3prime_exons = transcripts[0].utr3_exons[::-1]
    for i in range(1, len(transcripts)):
        t = transcripts[i]
        if t.strand == 1: exons = [[x[1],x[0]] for x in t.utr3_exons[::-1]]
        else: exons = t.utr3_exons[::-1]
        if most3prime_exons == exons: ret.append(t)
        elif (t.strand == 1 and exons > most3prime_exons) or (t.strand == -1 and exons < most3prime_exons):
            ret = [t]
            most3prime_exons = exons
    return ret
