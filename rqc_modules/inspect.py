def inspect(args):
    print("inspecting...")

    read_id = args.inputs[0]
    alignments = []

    for i in range(1, len(args.inputs), 2):
        label = args.inputs[i]
        filename = args.inputs[i+1]

        samfile = pysam.AlignmentFile(filename, 'rb')
        iter = samfile.fetch()

        for x in iter:
            if (x.query_name == read_id):
                alignments.append(x)

        samfile.close()

    for a in alignments:
        print(a)

