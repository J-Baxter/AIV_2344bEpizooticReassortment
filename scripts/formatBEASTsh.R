
# Set the file path
file_path <- "./scripts/eddie/beast_template.sh"
xmlfiles <- list.files('./data/alignments/subsampled_alignments/2024Jan24/xml/')

# Read the content of the file
lines <- readLines(file_path)

# Replace "XXX" with your new string


new_lines <- lapply(xmlfiles, function(x) gsub("XXX", x, lines))

# Write the modified content back to the file

shellnames <- paste0('./data/alignments/subsampled_alignments/2024Jan24/xml/', gsub('.xml$', '.sh', xmlfiles))

mapply(
  writeLines,
  new_lines,
  shellnames
)
