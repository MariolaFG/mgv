import sys
from ete2 import Tree, random_color, add_face_to_node, TreeStyle, TextFace, RectFace, AttrFace
import random
from collections import defaultdict
import argparse

import ete2
print ete2.__file__

def layout(node):
    node.img_style["size"] = 0
    node_name = node.name.replace(".bam", "")
    if node.is_leaf():
        col_country = header2column[sample_info["country"].get(node_name, "")]
        col_city = header2column[sample_info["city"].get(node_name, "")]
        ssize = 40
        countryF = RectFace(120, ssize, column_color[col_country], column_color[col_country])
        cityF = RectFace(120, ssize, column_color[col_city], column_color[col_city])

        if node_name in OBESITY:
            if OBESITY[node_name] == "obese":
                ocolor = "red"
            else:
                ocolor = "black"
            add_face_to_node(TextFace(OBESITY[node_name], fsize=18, fgcolor=ocolor), node, header2column["obesity"] , position="aligned")

        if node_name in DIABETES:
            ocolor = "purple"
            add_face_to_node(TextFace(DIABETES[node_name], fsize=18, fgcolor=ocolor), node, header2column["diabetes"] , position="aligned")

            
        add_face_to_node(countryF, node, col_country, position="aligned")
        add_face_to_node(cityF, node, col_city, position="aligned")
        add_face_to_node(TextFace("%s" %(node.name[:25]), fsize=20), node, 0, position='branch-right')
    else:
        if args.collapse_identical and len(node.children) > 50:
            node.img_style["draw_descendants"] = False
            node.img_style["size"] = 20
            node.img_style["shape"] = "square"
            node.img_style["hz_line_width"] = 5
            node.dist = 0.0001
            add_face_to_node(TextFace("%d samples without data" %len(node.children), fsize=20), node, 0, position='branch-right')
            
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("tree_files", metavar='tree_files', nargs='+')
    parser.add_argument('--show', dest='show', action='store_true')
    parser.add_argument('--collapse', dest='collapse_identical', action='store_true')
    args = parser.parse_args()
    print args
    
    sample_info = defaultdict(dict)
    for line in open('1086_sample_names_origin.txt'):
        name, city, country = line.strip().split('\t')
        sample_info["country"][name] = country
        sample_info["city"][name] = "%s - %s" %(country, city)

    OBESITY = {}
    DIABETES = {}
    for line in open('mh.meta'):
        fields = line.split('\t')
        OBESITY[fields[0].strip()]= fields[4].strip()
        DIABETES[fields[0].strip()]= fields[7].strip().replace("NA", "")
        
        
    column_color = ["black", "black"]
    column_header = ["obesity", "diabetes"]
    spaciators = set()
    for info_channel in sample_info.keys():
        unique_values = set(sample_info[info_channel].values())
        colors = random_color(num=len(unique_values))
        sorted_names = sorted(unique_values)
        column_header.extend(sorted_names)
        column_color.extend(colors)

        column_header.append("")
        column_color.append("black")
        spaciators.add(len(column_header)-1)
        
    header2column = dict([(name, i) for i, name in enumerate(column_header)])
        
    ts = TreeStyle()
    ts.mode = 'r'
    ts.draw_guiding_lines = False
    ts.show_leaf_name = False
    ts.force_topology = False
    ts.layout_fn = layout
    ts.tree_width = 800
    ts.draw_aligned_faces_as_table = True
    
    
    for i, name in enumerate(column_header):
        if name:
            headerF = TextFace(str(name), fgcolor=column_color[i], fsize=40)
            headerF.rotation = -85
        else:
            headerF = RectFace(300, 5, "white", "white")
        ts.aligned_header.add_face(headerF, i)

    #tree_files = sys.argv[1:]        
    for treefile in args.tree_files:
        output = treefile + '.png'
        print 'rendering', output
        try:
            t = Tree(open(treefile).read().replace('|', ','))
        except Exception, e:
            print e, treefile
        else:
            t.set_outgroup(t.get_midpoint_outgroup())
            t.sort_descendants()
            for n in t.traverse():
                if (n.children) > 2:
                    random.shuffle(n.children)                
                
                
            t.dist = 0.00001
            if args.show:
                t.show(tree_style=ts)
            else:
                t.render(output, tree_style=ts, w=1600)
                #t.render(treefile+'.pdf', tree_style=ts)
