# Edits a beauti XML file to remove clock and tree models to fit trait analyses
# to empirical trees.
#
# Arguments: 1) beauti XML that includes all your trait analysis
#            2) file path to posterior distribution of phylogenetic trees
#            (as obtained from BEAST)
#
# This script works reasonably well on typical BEAST parameterisations and model
# selection. If you have included advanced models in your Beauti XML, you may
# need to supplement 'elements.py'

# Copyright (c) 2024 James Baxter

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# import modules
import xml.etree.ElementTree as ET
import argparse
import re
import elements


def parse_args():
    parser = argparse.ArgumentParser(description="a script to do stuff")
    parser.add_argument("xml_path")
    parser.add_argument("empiricaltree_path")
    args = parser.parse_args()
    return args


def parse_query(key):
    q = ['', '', '']

    if key[1] and key[2] != '':
        q[1] = f"[@{key[1]}="

        if all(x != "" for x in key):
            q[2] = f"'{key[2]}']..."
        else:
            q[2] = f"'{key[2]}']"
    else:
        q[1] = key[1]
        q[2] = key[2]

    if key[0] == '':
        q[0] = '*'
    else:
        q[0] = key[0]

    query = f".//{q[0]}{q[1]}{q[2]}"
    return query


def remove_element(root, query):
    element_ids = root.findall(query)

    if isinstance(element_ids, list):
        for element in element_ids:
            for parent in root.iter():
                for child in list(parent):
                    if child is element:
                        parent.remove(child)

    else:
        for element in element_ids:
            root.remove(element)

    return root



#def remove_id(root, id_value):
   # element_ids = root.findall(f".//*[@id='{id_value}']")

   # for element in element_ids:
     #   root.remove(element)

   # return root


#def remove_parameters(root, id_ref):
    # element_ids = root.findall(f".//operators//parameter[@idref='{id_ref}']...")
    # could be generalised with f".//{element_name}[@{attrib_name}='{attrib_value}']"
    #element_ids = root.findall(f".//parameter[@idref='{id_ref}']...")

    #for element in element_ids:
      #  for parent in root.iter():
         #   for child in list(parent):
         #       if child is element:
          #          parent.remove(child)
   # return root


#def remove_others(root, element_name):
    #elements_named = root.findall(f".//{element_name}")

   # for element in elements_named:
       # for parent in root.iter():
          #  for child in list(parent):
              #  if child is element:
                 #   parent.remove(child)
    #return root


def add_empiricaltree(root, tree_path):

    new_sub1 = ET.Element('empiricalTreeDistributionModel')
    new_sub1.set("id", "treeModel")
    new_sub1.set("fileName", tree_path)
    new_sub1a = ET.SubElement(new_sub1, 'taxa')
    new_sub1a.set("idref", "taxa")
    position = root.find('taxa')
    root.insert(root.getchildren().index(position) + 1, new_sub1)

    new_sub2 = ET.Element('statistic')
    new_sub2.set("id", "treeModel.currentTree")
    new_sub2.set("name", 'Current Tree')
    new_sub2a = ET.SubElement(new_sub2, 'empiricalTreeDistributionModel')
    new_sub2a.set("idref", "treeModel")
    position = root.find('taxa')
    root.insert(root.getchildren().index(position) + 2, new_sub2)

    target_element = root.find(".//operators")
    new_sub3 = ET.SubElement(target_element, 'empiricalTreeDistributionOperator')
    new_sub3.set("weight", "3")
    new_sub3a = ET.SubElement(new_sub3, 'empiricalTreeDistributionModel')
    new_sub3a.set("idref", "treeModel")

    return root


def main():
    inputs = parse_args()
    xml = ET.parse(inputs.xml_path)
    root = xml.getroot()
    file_path = re.sub('.xml', '_empiricaltree.xml', inputs.xml_path)

    for key in keys:
        query = parse_query(key)
        root = remove_element(root, query)

    #for item in elements.id_to_remove:
     #   root = remove_id(root, item)

    #for item in elements.params_to_remove:
      #  root = remove_parameters(root, item)

    #for item in elements.others_to_remove:
        #root = remove_others(root, item)

    root = add_empiricaltree(root, inputs.tree_path)

    root.write(file_path)



if __name__ == "__main__":
    main()