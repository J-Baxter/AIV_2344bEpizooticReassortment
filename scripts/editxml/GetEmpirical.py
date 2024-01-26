#  import module
import xml.etree.ElementTree as ET
import argparse
from re import sub


# Parse beauty XML
def parse_args():
    parser=argparse.ArgumentParser(description="a script to do stuff")
    parser.add_argument("xml_path")
    parser.add_argument("empiricaltree_path")
    args=parser.parse_args()
    return args


def remove_id(xml, id_value):
    root = xml.getroot()
    elements_with_id = root.findall(f".//*[@id='{id_value}']")

    for element in elements_with_id:
        root.remove(element)

    return xml

test = remove_id(root, 'alignment')


def add_empiricaltree(xml, tree_path):
        # Create
        empiricalTreeDistributionModel
        {'id': 'treeModel', 'fileName': 'HA-H5_mper3_SRD06_relaxLn_skygrid6-72_2.combined.trees.txt'}
        # <empiricalTreeDistributionModel id="treeModel" fileName="HA-H5_mper3_SRD06_relaxLn_skygrid6-72_2.combined.trees.txt">
        #		<taxa idref="taxa"/>
        #	</empiricalTreeDistributionModel>
        #
        #	<statistic id="treeModel.currentTree" name="Current Tree">
        #		<empiricalTreeDistributionModel idref="treeModel"></empiricalTreeDistributionModel>
        #	</statistic>

def make_outputfile(xml_path):
    output_path = sub('.xml', '_empiricaltree.xml', xml_path)
    return output_path

def main():
    inputs = parse_args()
    xml = ET.parse(inputs.xml_path)
    file_path = make_outputfile(inputs.xml_path)

    for element_id in element_id_list:
        remove_id(xml, element_id)

    add_empiricaltree(xml, inputs.tree_path)


    xml.write(file_path)









# Remove attributes associated with alignment, clock and tree models

# Absolute removals:
alignment {'id': 'alignment', 'dataType': 'nucleotide'}
#mergePatterns {'id': 'CP1+2.patterns'}
#patterns {'id': 'CP3.patterns', 'from': '3', 'every': '3', 'strip': 'false'}
#constantSize {'id': 'initialDemo', 'units': 'years'}
#coalescentSimulator {'id': 'startingTree'}
#treeModel {'id': 'treeModel'}
#treeLengthStatistic {'id': 'treeLength'}
#tmrcaStatistic {'id': 'age(root)', 'absolute': 'true'}
#gmrfSkyGridLikelihood {'id': 'skygrid'}
#discretizedBranchRates {'id': 'branchRates'}
#rateCovarianceStatistic {'id': 'covariance', 'name': 'covariance'}
#HKYModel {'id': 'CP1+2.hky'}
#HKYModel {'id': 'CP3.hky'}
#compoundParameter {'id': 'allMus'}
#treeDataLikelihood {'id': 'treeLikelihood', 'useAmbiguities': 'false'}

# Specific ID's
#siteModel {'id': 'CP1+2.siteModel'} {'id': 'CP3.siteModel'}
#rateStatistic {'id': 'meanRate', 'name': 'meanRate', 'mode': 'mean', 'internal': 'true', 'external': 'true'}
#rateStatistic {'id': 'coefficientOfVariation', 'name': 'coefficientOfVariation', 'mode': 'coefficientOfVariation', 'internal': 'true', 'external': 'true'}




# Keep
taxa {'id': 'taxa'}
generalDataType {'id': 'Subtype.dataType'}
generalDataType {'id': 'Host.dataType'}
attributePatterns {'id': 'Subtype.pattern', 'attribute': 'Subtype'}
attributePatterns {'id': 'Host.pattern', 'attribute': 'Host'}
empiricalTreeDistributionModel {'id': 'treeModel', 'fileName': 'HA-H5_mper3_SRD06_relaxLn_skygrid6-72_2.combined.trees.txt'}
statistic {'id': 'treeModel.currentTree', 'name': 'Current Tree'}
strictClockBranchRates {'id': 'Subtype.branchRates'}
rateStatistic {'id': 'Subtype.meanRate', 'name': 'Subtype.meanRate', 'mode': 'mean', 'internal': 'true', 'external': 'true'}
strictClockBranchRates {'id': 'Host.branchRates'}
rateStatistic {'id': 'Host.meanRate', 'name': 'Host.meanRate', 'mode': 'mean', 'internal': 'true', 'external': 'true'}
multivariateDiffusionModel {'id': 'latlon.diffusionModel'}
multivariateWishartPrior {'id': 'latlon.precisionPrior', 'df': '2'}
multivariateTraitLikelihood {'id': 'latlon.traitLikelihood', 'traitName': 'latlon', 'useTreeLength': 'true', 'scaleByTime': 'true', 'reportAsMultivariate': 'true', 'reciprocalRates': 'true', 'integrateInternalTraits': 'true'}
correlation {'id': 'latlon.correlation', 'dimension1': '1', 'dimension2': '2'}
matrixInverse {'id': 'latlon.varCovar'}
continuousDiffusionStatistic {'id': 'latlon.diffusionRate', 'greatCircleDistance': 'true'}
generalSubstitutionModel {'id': 'Subtype.model'}
sumStatistic {'id': 'Subtype.nonZeroRates', 'elementwise': 'true'}
productStatistic {'id': 'Subtype.actualRates', 'elementwise': 'false'}
siteModel {'id': 'Subtype.siteModel'}
generalSubstitutionModel {'id': 'Host.model', 'randomizeIndicator': 'false'}
sumStatistic {'id': 'Host.nonZeroRates', 'elementwise': 'true'}
productStatistic {'id': 'Host.actualRates', 'elementwise': 'false'}
siteModel {'id': 'Host.siteModel'}
ancestralTreeLikelihood {'id': 'Subtype.treeLikelihood', 'stateTagName': 'Subtype.states', 'useUniformization': 'true', 'saveCompleteHistory': 'false', 'logCompleteHistory': 'false'}
ancestralTreeLikelihood {'id': 'Host.treeLikelihood', 'stateTagName': 'Host.states', 'useUniformization': 'true', 'saveCompleteHistory': 'false', 'logCompleteHistory': 'false'}
operators {'id': 'operators', 'optimizationSchedule': 'log'}
mcmc {'id': 'mcmc', 'chainLength': '1100000', 'autoOptimize': 'true', 'operatorAnalysis': 'HA-H5_mper3_empirical_Subtype_Host_latlon_0425.ops.txt'}
report {}


newtree = ET.parse('./scripts/eddie/editxml/HA'
newroot = tree.getroot()

id_value = ['alignment', 'dataType', 'nucleotide', 'CP1+2.patterns', ]

target_elements = newroot.findall(f".//*[@id='{id_value}']")

for element in target_elements:
    print(ET.tostring(element, encoding='utf-8').decode('utf-8'))


def remove_attribute(root, element_name):
    for target in root.findall(element_name):
        parent = target.getparent()
        parent.remove(element)

    return parent



for target in root.findall('alignment'):
    parent = target.getparent()






if __name__ == "__main__":
    main()

    for neighbor in tree.iter('neighbor'):
        print(neighbor.attrib)