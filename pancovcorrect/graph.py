"""
Function use to convert assembly in graph
"""


def asm2gfa(input_path, output_path, k):
    """
    Convert bcalm assemblie in gfa format
    """

    with open(input_path) as input_file:
        # name stores the id of the unitig
        # optional is a list which stores all the optional tags of a segment
        # links stores all the link information about a segment
        name = ""
        optional = []
        links = []
        graph = open(output_path, "w")
        # adding Header to the file
        graph.write("H\tVN:Z:1.0\tks:i:%d\n" % k)  # includes the k-mer size

        # firstLine is for implemetation purpose so that we don't add some
        # garbage value to the output file.
        first_line = 0
        # segment stores the segment till present, in a fasta file,
        # segment can be on many lines, hence we need to get the whole segment
        # from all the lines
        segment = ""
        for line in input_file:
            line = line.replace("\n", "")
            if line[0] != ">":
                # segment might be in more than one line, hence we get the
                # whole segment first, and then put it in the GFA file.
                segment += line
            if line[0] == ">":
                if (
                    first_line != 0
                ):  # if it's not the firstline in the input file, we store
                    # the input in GFA format in the output file
                    __write_segment(name, segment, optional, graph, links)
                    segment = ""

                first_line = 1
                # once the previous segment and it's information has been
                # stored, we start the next segment and it's information
                line = line.split(" ")
                name = line[0][1:]  # get the id
                optional = []
                links = []
                # we skip the first value because the first value is ">ID"
                for i in range(1, len(line)):
                    # we need this because the line can end with a space, hence
                    # we get one extra value in our list.
                    if line[i] == "":
                        continue
                    if (
                        line[i][0:2] == "MA"
                    ):  # previous bcalm2 versions had "MA=[xxx]" optional tag
                        # as well, kept it just for compatibility, and reformated
                        optional.append(line[i][0:2] + ":f:" + line[i][2:])
                    elif line[i][0:2] == "L:":  # for links
                        link = line[i].split(":")
                        links.append(
                            "L\t"
                            + name
                            + "\t"
                            + link[1]
                            + "\t"
                            + link[2]
                            + "\t"
                            + link[3]
                            + "\t"
                            + str(k - 1)
                            + "M\n"
                        )
                    else:  # all the other optional tags
                        optional.append(line[i])

        # we will miss the last one, because it won't go into the if condition
        # - if(line[0]==">") and hence won't add the segment to the file.
        __write_segment(name, segment, optional, graph, links)

        graph.close()


def __write_segment(name, segment, optional, graph, links):
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
    graph.write(add.strip() + "\n")
    for (
        j
    ) in links:  # adding all the links of the current segment to the GFA file
        graph.write(j)
