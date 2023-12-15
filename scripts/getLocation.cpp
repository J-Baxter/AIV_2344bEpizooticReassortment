#include <Rcpp.h>
#include <unordered_map>

// Define a type for the atlas information
struct Atlas {
  std::string continent;
  std::string country;
  std::string subdivision;
  //latlong
};

// Function to create a hash table of countries
std::unordered_map<std::string, Atlas> createAtlasHashTable() {
  std::unordered_map<std::string, Atlas> atlasTable;
  
  // Add atlas information
  atlasTable["common ostrich"]= {"Struthioniformes", "Struthionidae", "struthio camelus"};
  atlasTable["somali ostrich"]= {"Struthioniformes", "Struthionidae", "struthio molybdophanes"};
  atlasTable["common/somali ostrich"]= {"Struthioniformes", "Struthionidae", "struthio camelus/molybdophanes"};
  atlasTable["greater rhea"]= {"Rheiformes", "Rheidae", "rhea americana"};
  atlasTable["lesser rhea"]= {"Rheiformes", "Rheidae", "rhea pennata"};
  atlasTable["tawny-breasted tinamou"]= {"Tinamiformes", "Tinamidae", "nothocercus julius"};
  atlasTable["highland tinamou"]= {"Tinamiformes", "Tinamidae", "nothocercus bonapartei"};
  atlasTable["hooded tinamou"]= {"Tinamiformes", "Tinamidae", "nothocercus nigrocapillus"};
  
  
  return atlasTable;
}

// Function to get atlas for a given location
Rcpp::List getAtlas(std::unordered_map<std::string, Atlas>& atlasTable, const std::string& commonName) {
  
  // Check if the location exists in the hash table
  auto it = atlasTable.find(commonName);
  if (it != atlasTable.end()) {
    Atlas queryInfo = it->second;
    
    // Create a list with atlas information
    Rcpp::DataFrame result = 
      Rcpp::DataFrame::create(Rcpp::Named("continent") = queryInfo.continent,
                              Rcpp::Named("country") = queryInfo.country,
                              Rcpp::Named("subdivision") = queryInfo.subdivision);//latlonh
    return result;
  } else {
    Rcpp::stop("Location not found in the hash table.");
  }
}

// [[Rcpp::export]]
Rcpp::List getLocation(const std::string& commonName) {
  // Create a hash table of bird atlas
  std::unordered_map<std::string, Atlas> atlasTable = createAtlasHashTable();
  
  // Get bird atlas for the specified common name
  return getAtlas(atlasTable, commonName);
}