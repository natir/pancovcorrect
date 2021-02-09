"""
Function use to convert assembly in graph
"""


def asm2gfa(input_path, output_path, k):
    """
    Convert bcalm assemblie in gfa format
    """

    with open(input_path) as f:
        # name stores the id of the unitig
        # optional is a list which stores all the optional tags of a segment
        # links stores all the link information about a segment
        name = ""
        optional = []
        links = []
        g = open(output_path, "w")
        # adding Header to the file
        g.write("H\tVN:Z:1.0\tks:i:%d\n" % k)  # includes the k-mer size

        # firstLine is for implemetation purpose so that we don't add some
        # garbage value to the output file.
        firstLine = 0
        # segment stores the segment till present, in a fasta file,
        # segment can be on many lines, hence we need to get the whole segment
        # from all the lines
        segment = ""
        for line in f:
            line = line.replace("\n", "")
            if line[0] != ">":
                # segment might be in more than one line, hence we get the
                # whole segment first, and then put it in the GFA file.
                segment += line
            if line[0] == ">":
                if (
                    firstLine != 0
                ):  # if it's not the firstline in the input file, we store
                    # the input in GFA format in the output file
                    __write_segment(name, segment, optional, g, links)
                    segment = ""

                firstLine = 1
                # once the previous segment and it's information has been
                # stored, we start the next segment and it's information
                a = line.split(" ")
                name = a[0][1:]  # get the id
                optional = []
                links = []
                # we skip the first value because the first value is ">ID"
                for i in range(1, len(a)):
                    # we need this because the line can end with a space, hence
                    # we get one extra value in our list.
                    if a[i] == "":
                        continue
                    if (
                        a[i][0:2] == "MA"
                    ):  # previous bcalm2 versions had "MA=[xxx]" optional tag
                        # as well, kept it just for compatibility, and reformated
                        optional.append(a[i][0:2] + ":f:" + a[i][2:])
                    elif a[i][0:2] == "L:":  # for links
                        b = a[i].split(":")
                        k1 = k - 1
                        links.append(
                            "L\t"
                            + name
                            + "\t"
                            + b[1]
                            + "\t"
                            + b[2]
                            + "\t"
                            + b[3]
                            + "\t"
                            + str(k1)
                            + "M\n"
                        )
                    else:  # all the other optional tags
                        optional.append(a[i])

        # we will miss the last one, because it won't go into the if condition
        # - if(line[0]==">") and hence won't add the segment to the file.
        __write_segment(name, segment, optional, g, links)

        g.close()


def __write_segment(name, segment, optional, g, links):
    add = ""
    add += "S\t"  # for segment
    add += name  # id of segment
    add += "\t"
    add += segment  # segment itself
    add += "\t"
    for i in optional:  # optional tags
        add += i
        add += "\t"
    # adding Segment to the file
    g.write(add.strip() + "\n")
    for (
        j
    ) in links:  # adding all the links of the current segment to the GFA file
        g.write(j)
