# Edits a beauti XML file to remove clock and tree models to fit trait analyses
# to empirical trees.
#
# Arguments: 1) beauti XML that includes all your trait analysis
#            2) file path to posterior distribution of phylogenetic trees
#            (as obtained from BEAST)
#
# This script works reasonably well using default BEAST parameterisations for clock and tree model.
# If you have included advanced models in your Beauti XML, you may need to supplement the keys list.

# NB: requires python 3.9 or later

# Copyright (c) 2024 James Baxter under GNU GENERAL PUBLIC LICENSE Version 3

# import modules
import xml.etree.ElementTree as ET
import argparse
import re

keys = [['alignment', '', ''],
        ['patterns', '', ''],
        ['constantSize', '', ''],
        ['coalescentSimulator', '', ''],
        ['treeModel', '', ''],
        ['treeLengthStatistic', '', ''],
        ['tmrcaStatistic', '', ''],
        ['coalescentLikelihood', '', ''],
        ['HKYModel', '', ''],
        ['statistic', '', ''],
        ['treeDataLikelihood', '', ''],
        ['upDownOperator', '', ''],
        ['subtreeSlide', '', ''],
        ['narrowExchange', '', ''],
        ['wideExchange', '', ''],
        ['wilsonBalding', '', ''],
        ['uniformOperator', '', ''],
        ['logNormalPrior', '', ''],
        ['dirichletPrior', '', ''],
        ['', 'id', 'default.branchRates'],
        ['', 'id', 'default.meanRate'],
        ['', 'id', 'mu'],
        ['', 'id', 'siteModel'],
        ['', 'idref', 'kappa'],
        ['', 'idref', 'frequencies'],
        ['', 'idref', 'default.clock.rate'],
        ['', 'idref', 'default.meanRate'],
        ['', 'idref', 'treeModel.allInternalNodeHeights'],
        ['', 'idref', 'treeModel.rootHeight'],
        ['', 'idref', 'constant.popSize'],
        ['', 'idref', 'default.branchRates'],
        ['', 'label', 'default.clock.rate'],
        ['', 'label', 'age(root)'],
        ['', 'tag', 'default.rate'],
        ['oneOnXPrior', '', '']]

def parse_args():
    parser = argparse.ArgumentParser(description="a script to do stuff")
    parser.add_argument("xml_path")
    parser.add_argument("tree_path")
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


def add_empiricaltree(root, tree_path):
    new_sub1 = ET.Element('empiricalTreeDistributionModel')
    new_sub1.set("id", "treeModel")
    new_sub1.set("fileName", tree_path)
    new_sub1a = ET.SubElement(new_sub1, 'taxa')
    new_sub1a.set("idref", "taxa")

    # after taxa but before 'strictClockBranchRates' = additional flexibility required
    position = root.find('taxa')
    root.insert(list(root).index(position) + 1, new_sub1)

    new_sub2 = ET.Element('statistic')
    new_sub2.set("id", "treeModel.currentTree")
    new_sub2.set("name", 'Current Tree')
    new_sub2a = ET.SubElement(new_sub2, 'empiricalTreeDistributionModel')
    new_sub2a.set("idref", "treeModel")
    position = root.find('taxa')
    root.insert(list(root).index(position) + 2, new_sub2)

    target_element = root.find(".//operators")
    new_sub3 = ET.SubElement(target_element, 'empiricalTreeDistributionOperator')
    new_sub3.set("weight", "3")
    new_sub3a = ET.SubElement(new_sub3, 'empiricalTreeDistributionModel')
    new_sub3a.set("idref", "treeModel")

    target_element = root.find(".//logTree")
    new_sub4 = ET.Element('treeModel')
    new_sub4.set("idref", "treeModel")
    target_element.insert(0, new_sub4)

    target_element = root.findall(".//ctmcScalePrior")
    if target_element:
        for target in target_element:
            target.insert(1, new_sub4)

    target_element = root.find(".//multivariateTraitLikelihood[multivariateDiffusionModel]")
    # Check if the element is found
    if target_element is not None:
        target_element.insert(1, new_sub4)

    target_element = root.findall(".//ancestralTreeLikelihood[attributePatterns]")
    # Check if the element is found
    if target_element:
        for target in target_element:
            target.insert(1, new_sub4)

    for rate_statistic in root.findall(".//rateStatistic[@mode='mean']"):
        rate_statistic.insert(0, new_sub4)

    return root


def main():
    inputs = parse_args()
    xml = ET.parse(inputs.xml_path)
    root = xml.getroot()
    file_path = re.sub('.xml', '_empiricaltree.xml', inputs.xml_path)

    for key in keys:
        query = parse_query(key)
        root = remove_element(root, query)

    operators_element = root.find(".//operators")

    # If <operators> element exists
    if operators_element:
        # Find all <scaleOperator> elements under <operators>
        scale_operators = operators_element.findall("./scaleOperator")

        # Iterate over each <scaleOperator> element
        for scale_operator in scale_operators:
            # Check if the <scaleOperator> element is empty (has no subelements)
            if len(scale_operator) == 0:
                # Remove the empty <scaleOperator> element from <operators>
                operators_element.remove(scale_operator)

        scale_operators = operators_element.findall("./deltaExchange")
        # Iterate over each <scaleOperator> element
        for scale_operator in scale_operators:
            # Check if the <scaleOperator> element is empty (has no subelements)
            if len(scale_operator) == 0:
                # Remove the empty <scaleOperator> element from <operators>
                operators_element.remove(scale_operator)

    mcmc_element = root.find(".//mcmc")

    if mcmc_element:
        joint_element = mcmc_element.find("./joint")
        if joint_element:
            # Find the <prior> element under <joint>
            prior_element = joint_element.find("./prior")
            # If <prior> element exists
            if prior_element:
                # Find the <ctmcScalePrior> element under <prior>
                ctmc_scale_prior_element = prior_element.find("./ctmcScalePrior")
                # If <ctmcScalePrior> element exists
                if ctmc_scale_prior_element:
                    # Find the <ctmcScale> element under <ctmcScalePrior>
                    ctmc_scale_element = ctmc_scale_prior_element.find("./ctmcScale")
                    # If <ctmcScale> element exists
                    if len(ctmc_scale_element) == 0:
                        # Remove <ctmcScale> from its parent <ctmcScalePrior>
                        ctmc_scale_prior_element.remove(ctmc_scale_element)

                if len(ctmc_scale_prior_element) == 0:
                # Remove <ctmcScale> from its parent <ctmcScalePrior>
                    prior_element.remove(ctmc_scale_prior_element)


    root = add_empiricaltree(root, inputs.tree_path)

    tree = ET.ElementTree(root)
    ET.indent(tree, space="\t", level=0)
    tree.write(file_path, encoding="utf-8")



if __name__ == "__main__":
    main()
