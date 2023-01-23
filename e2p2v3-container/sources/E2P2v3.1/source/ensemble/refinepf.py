import re
# from argparse import ArgumentParser


def read_pf(fp):
    id, attr = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith('ID\t'):
            if id: yield (id, '\n'.join(attr))
            id, attr = line, []
        else:
            attr.append(line)
    if id: yield (id, '\n'.join(attr))


def remove_empty_from_pf(pf_path):
    refined_dict = {}
    exist_empty = False
    empty_count = 0
    with open(pf_path, 'r') as fp:
        print('Opening pf file:\t' + pf_path)
        for e2p2_id, e2p2_attr in read_pf(fp):
            metacyc_attr = re.findall(r'METACYC\t[^\s]+\n', e2p2_attr)
            if len(metacyc_attr) > 0:
                refined_dict.setdefault(e2p2_id, e2p2_attr)
            else:
                # print("Empty Entry: %s", e2p2_id.replace('ID\t', '').strip())
                empty_count += 1
                exist_empty = True
    if exist_empty:
        print('Removed Empty IDs:\t' + str(empty_count))
        with open(pf_path, 'w') as fp:
            print('Rewriting pf file:\t' + pf_path)
            for e2p2_id in sorted(refined_dict.keys()):
                fp.write(e2p2_id + '\n' + refined_dict[e2p2_id] + '\n')


if __name__ == '__main__':
    pass
    # parser = ArgumentParser()
    # parser.add_argument('input_path')

    # args = parser.parse_args()
    # remove_empty_from_pf(args.input_path)

