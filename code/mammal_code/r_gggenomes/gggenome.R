
library(gggenomes)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(ggpattern) 

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("rna-directed rna polymerase", "RdRp", interpro_description, ignore.case = TRUE) %>%
           gsub("rna dependent rna polymerase", "RdRp", ., ignore.case = TRUE) %>%
           gsub("rna-dependent rna polymerase", "RdRp", ., ignore.case = TRUE) %>%
           gsub("RNA-directed RNA polymerase L, C-terminal", "RdRp", ., ignore.case = TRUE) %>%
           gsub("dna/rna polymerase superfamily", "RdRp", ., ignore.case = TRUE) %>%
           gsub("rdrp,?\\s*c?-?terminal\\s*(,\\s*virus)?", "RdRp", ., ignore.case = TRUE) %>%
           gsub("rdrp", "RdRp", ., ignore.case = TRUE) %>%
           gsub("RdRp,\\s*catalytic domain", "RdRp", ., ignore.case = TRUE) %>%
           # Manually replace new patterns
           gsub("RdRp beta-chain", "RdRp", ., ignore.case = TRUE) %>%
           gsub("Mononegavirales RdRp catalytic domain", "RdRp, Mononegavirales", ., ignore.case = TRUE) %>%
           gsub("RdRp L", "RdRp", ., ignore.case = TRUE) %>%
           gsub("RdRpdomain", "RdRp", ., ignore.case = TRUE) %>%
           gsub("RdRp,\\s*virus", "RdRp", ., ignore.case = TRUE))


selected_rows <- selected_rows %>%
  mutate(interpro_description = ifelse(grepl("helicase|dead|p-loop", interpro_description, ignore.case = TRUE), 
                                       "Helicase", 
                                       interpro_description))

selected_rows$interpro_description <- gsub("(?i)macro domain(-like)?", "MD", selected_rows$interpro_description, perl = TRUE)

library(dplyr)

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub(
    "^Nucleocapsid protein, (N-terminal|C-terminal)", "N-protein",
    interpro_description, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Nucleocapsid \\(N\\) protein, (N-terminal domain|C-terminal domain)(,\\s*(.*))?$",
    "N-protein \\3",
    interpro_description, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Nucleocapsid protein, (coronavirus)(,\\s*(.*))?$",
    "N-protein \\1",
    interpro_description, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Nucleocapsid protein, (.*),\\s*(.*)$",
    "N-protein \\2",
    interpro_description, perl = TRUE
  ))



selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub(
    "^Arterivirus papain-like cysteine protease alpha \\(PCPalpha\\) domain( superfamily)?$", "PCPalpha",
    interpro_description, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Arterivirus papain-like cysteine protease beta \\(PCPbeta\\) domain$", "PCPbeta",
    interpro_description, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Arterivirus papain-like cysteine protease beta \\(PCPbeta\\) as PCPalpha and Beta$", "PCPalpha and Beta",
    interpro_description, perl = TRUE
  ))




selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub(
    "^Pestivirus envelope glycoprotein E2(, domain (A|B|D))?$", "E2",
    interpro_description, perl = TRUE
  ))


selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub(
    "^Spike glycoprotein S1, coronavirus$", "Spike1 coronavirus",
    interpro_description, ignore.case = TRUE, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Spike glycoprotein S2, coronavirus(, heptad repeat (1|2))?(, C-terminal)?( superfamily)?$", "Spike2 coronavirus",
    interpro_description, ignore.case = TRUE, perl = TRUE
  ))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub(
    "^M matrix/glycoprotein, (alphacoronavirus|coronavirus)$", "MP \\1",
    interpro_description, ignore.case = TRUE, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Pneumovirus matrix 2-1$", "MP Pneumovirus",
    interpro_description, ignore.case = TRUE, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Vesiculovirus matrix$", "MP Vesiculovirus",
    interpro_description, ignore.case = TRUE, perl = TRUE
  )) %>%
  mutate(interpro_description = gsub(
    "^Viral matrix protein(, (C-terminal domain|N-terminal domain))?$", "MP",
    interpro_description, ignore.case = TRUE, perl = TRUE
  ))

selected_rows$interpro_description <- gsub("Pneumovirus M2-2", "MP pneumovirus", selected_rows$interpro_description)


# Standardize the descriptions
selected_rows$interpro_description <- gsub(
  "^Zinc finger, CCCH-type.*$", 
  "Zinc finger, CCCH-type",
  selected_rows$interpro_description,
  ignore.case = TRUE
)

# Group matrix descriptions
selected_rows$interpro_description <- gsub(
  "^matrix\\s*(alphacoronavirus|coronavirus)?$",
  "MP coronavirus",
  selected_rows$interpro_description,
  ignore.case = TRUE
)

# Replace 'Non-structural protein' with 'NSP'
selected_rows$interpro_description <- gsub(
  "^Non-structural protein\\s+", 
  "NSP", 
  selected_rows$interpro_description,
  ignore.case = TRUE
)

# Ensure there is no double 'NSP' in the description
selected_rows$interpro_description <- gsub(
  "(NSP)\\1", 
  "\\1", 
  selected_rows$interpro_description
)


selected_rows$interpro_description <- gsub("alphacoronavirus", "coronavirus", selected_rows$interpro_description, ignore.case = TRUE)

# Replace "RNA polymerase, N-terminal, coronavirus" with "RdRp, coronavirus" in the interpro_description column
selected_rows$interpro_description <- gsub("RNA polymerase, N-terminal, coronavirus", "RdRp, coronavirus", selected_rows$interpro_description, fixed = TRUE)
# Replace the specified phrases with "NSP10, coronavirus"
selected_rows$interpro_description <- gsub("RNA synthesis protein NSP10 superfamily, coronavirus", "NSP10, coronavirus", selected_rows$interpro_description, fixed = TRUE)
selected_rows$interpro_description <- gsub("RNA synthesis protein NSP10, coronavirus", "NSP10, coronavirus", selected_rows$interpro_description, fixed = TRUE)
selected_rows$interpro_description <- gsub("Nonstructural protein (\\d+).*?,\\s*coronavirus.*", "NSP\\1, coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Nonstructural protein (\\d+).*?,\\s*coronavirus-like", "NSP\\1, coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Nonstructural protein (\\d+).*?,\\s*alpha/beta-coronavirus", "NSP\\1, coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Nonstructural protein (\\d+).*?,\\s*alpha/betacoronavirus", "NSP\\1, coronavirus", selected_rows$interpro_description)

# Clean the descriptions
selected_rows$interpro_description <- gsub("NSP(\\d+).*?,\\s*coronavirus.*", "NSP\\1, coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("NSP(\\d+).*?,\\s*coronavirus-like", "NSP\\1, coronavirus", selected_rows$interpro_description)

# Clean the descriptions
selected_rows$interpro_description <- gsub("Peptidase C(\\d+).*?,\\s*coronavirus", "Pep C\\1, coronavirus", selected_rows$interpro_description)


# Apply specific replacements
selected_rows$interpro_description <- gsub("Arterivirus Nsp11 N-terminal/Coronavirus NSP15 middle domain", "NSP15, coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Coronavirus replicase NSP15, N-terminal oligomerization", "NSP15, coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Coronavirus replicase NSP3, C-terminal", "NSP3, coronavirus", selected_rows$interpro_description)

selected_rows$interpro_description <- gsub("Spike glycoprotein S2 superfamily, coronavirus", "Spike2 coronavirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Coronavirus spike glycoprotein S1, C-terminal", "Spike1 coronavirus", selected_rows$interpro_description)


selected_rows$interpro_description <- gsub("NSP15, middle domain superfamily", "NSP15", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("NSP11, NendoU domain, arterivirus", "NSP11, arterivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("NSP11, N-terminal, arterivirus", "NSP11, arterivirus", selected_rows$interpro_description)
# Apply specific replacements
selected_rows$interpro_description <- gsub("Arterivirus polyprotein, nsp2, immunogenic region", "NSP2, arterivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Arterivirus nsp7 alpha superfamily", "NSP7, arterivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Arterivirus NSP4 peptidase domain", "NSP4, arterivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Arterivirus Nsp2, peptidase C33", "NSP2, arterivirus", selected_rows$interpro_description)

selected_rows$interpro_description <- gsub("NSP3, X-domain-like", "NSP3", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Nonstructural protein 2, HCoV-229E-like", "NSP2", selected_rows$interpro_description)

selected_rows$interpro_description <- gsub("Pestivirus nonstructural protein 2", "NSP2, pestivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Nonstructural protein 10, 1B domain, arterivirus", "NSP10, arterivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Nonstructural protein 10, zinc-binding domain, arterivirus", "NSP10, arterivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Arterivirus nonstructural protein 7 alpha", "NSP7, arterivirus", selected_rows$interpro_description)

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Arterivirus GP2a envelope protein" = "EnvGP2a",
                                                  "Arterivirus GP3 envelope glycoprotein" = "EnvGP3",
                                                  "Arterivirus GP5 envelope glycoprotein" = "EnvGP5")))

selected_rows$interpro_description <- gsub("Peptidase C53, pestivirus Npro", "Pep C53, pestivirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Peptidase C53, pestivirus Npro, interaction domain", "Pep C53, pestivirus", selected_rows$interpro_description)

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Peptidase S1, PA clan, chymotrypsin-like fold" = "Pep S1",
                                                  "Peptidase S1, PA clan" = "Pep S1")))
selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Pestivirus NS2, peptidase C74" = "Pep C74, pestivirus",
                                                  "Pestivirus NS3, peptidase S31" = "Pep S31, pestivirus")))

selected_rows$interpro_description <- selected_rows$interpro_description %>%
  gsub("3a-like viroporin, cytosolic domain, alpha/betacoronavirus", "3a-like VP, coronavirus", .) %>%
  gsub("3a-like viroporin, transmembrane domain, alpha/betac", "3a-like VP, coronavirus", .)

# Perform the replacement
selected_rows$interpro_description <- selected_rows$interpro_description %>%
  gsub("RdRp, C-terminal", "RdRp", .)

selected_rows$interpro_description <- selected_rows$interpro_description %>%
  gsub("RdRp, C-terminal", "RdRp", .)

# Assuming your data is in a vector called selected_rows$interpro_description
selected_rows$interpro_description <- gsub("Ribonuclease T2-like superfamily", "RNase T2", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Ribonuclease T2, His active site 2", "RNase T2", selected_rows$interpro_description)

# Replace specific terms in the interpro_description column
selected_rows$interpro_description <- gsub("Nidovirus 3'-5' exoribonuclease domain", "Exo RNase, nidovirus", selected_rows$interpro_description)
selected_rows$interpro_description <- gsub("Endoribonuclease EndoU-like", "Endo RNase", selected_rows$interpro_description)

# Replace specific term in the interpro_description column
selected_rows$interpro_description <- gsub("Equine arteritis virus peptidase S32", "Pep S32, Equine arteritis", selected_rows$interpro_description)
# Replace specific term in the interpro_description column
selected_rows$interpro_description <- gsub("Pep C53, pestivirus, interaction domain", "Pep C53, pestivirus", selected_rows$interpro_description)
# Replace specific term in the interpro_description column
selected_rows$interpro_description <- gsub("Hepatitis E virus, cysteine peptidase", "Pep, Hep E virus", selected_rows$interpro_description)
# Replace specific term in the interpro_description column
selected_rows$interpro_description <- gsub("3a-like VP, coronavirusoronavirus", "3a-like VP, coronavirus", selected_rows$interpro_description)

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace(interpro_description, 
                                            "Paramyxovirus nucleocapsid protein", 
                                            "N-protein"))


selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace(interpro_description, 
                                            "Nidovirus RdRp-associated nucleotidyl transferase (NiRAN) domain", 
                                            "NiRAN"))

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Alphavirus-like methyltransferase \\(MT\\) domain" = "MT",
                                                  "RdRp methyltransferase domain, rhabdovirus" = "RdRp",
                                                  "RdRp methyltransferase domain, paramyxovirus" = "RdRp",
                                                  "Nidovirus 2-O-methyltransferase" = "MT",
                                                  "Mononegavirus L protein 2-O-ribose methyltransferase" = "MT",
                                                  "Nidovirus RdRp-associated nucleotidyl transferase \\(NiRAN\\) domain" = "NiRAN")))

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_trim(interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Nidovirus\\s+RdRp-associated\\s+nucleotidyl\\s+transferase\\s*\\(NiRAN\\)\\s*domain" = "NiRAN")))

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Rhabdovirus\\s+glycoprotein" = "GP",
                                                  "Major\\s+surface\\s+glycoprotein\\s+G" = "GP",
                                                  "Precursor\\s+fusion\\s+glycoprotein\\s+F0,\\s+Paramyxoviridae" = "GP")))

selected_rows <- selected_rows %>%
  mutate(interpro_description = str_replace_all(interpro_description, 
                                                c("Pep\\s+S1" = "Pep",
                                                  "Pep\\s+C30,\\s+coronavirus" = "Pep",
                                                  "Pep\\s+C16,\\s+coronavirus" = "Pep",
                                                  "Pep\\s+S31,\\s+pestivirus" = "Pep",
                                                  "Pep\\s+C53,\\s+pestivirus" = "Pep",
                                                  "Pep\\s+C74,\\s+pestivirus" = "Pep",
                                                  "Pep,\\s+Hep\\s+E\\s+virus" = "Pep",
                                                  "Pep\\s+S32,\\s+Equine\\s+arteritis" = "Pep",
                                                  "Replicase\\s+polyprotein\\s+1ab,\\s+peptidase\\s+C33-associated" = "Pep",
                                                  "Pep\\s+domain" = "Pep")))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("S-adenosyl-L-methionine-dependent methyltransferase superfamily", "MT", interpro_description))


selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Arterivirus papain-like cysteine protease beta \\(PCPbeta\\) domain superfamily", "PCP(a/ß)", interpro_description)) %>%
  mutate(interpro_description = gsub("PCPalpha", "PCP(a/ß)", interpro_description)) %>%
  mutate(interpro_description = gsub("PCPbeta", "PCP(a/ß)", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("RNase T2", "RNase", interpro_description)) %>%
  mutate(interpro_description = gsub("Endo RNase", "RNase", interpro_description)) %>%
  mutate(interpro_description = gsub("Exo RNase, nidovirus", "RNase", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Viral coat protein subunit", "coat protein", interpro_description)) %>%
  mutate(interpro_description = gsub("Picornavirus/Calicivirus coat protein", "coat protein", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("MP pneumovirus", "MP", interpro_description)) %>%
  mutate(interpro_description = gsub("MP Vesiculovirus", "MP", interpro_description)) %>%
  mutate(interpro_description = gsub("MP coronavirus", "MP", interpro_description)) %>%
  mutate(interpro_description = gsub("MP Pneumovirus", "MP", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("RdRp,?\\s*(bacteriophage|hepatitis C virus|Mononegavirales|orbiviral|coronavirus|luteovirus)?", "RdRp", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = ifelse(
    grepl("^NSP\\d+", interpro_description), 
    gsub(",.*", "", interpro_description), 
    interpro_description
  ))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Spike1 coronavirus|Spike2 coronavirus", "Spike(1,2)", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Haemagglutinin/haemagglutinin-neuraminidase, paramyxovirus", "HN", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Papain-like protease, thumb domain superfamily, coronavirus", "PLpro", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Serine protease, chymotrypsin-like serine protease, C-terminal", "CT-SP", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("N-protein coronavirus", "N-protein", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Pectin lyase fold/virulence factor", "PL", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Reverse transcriptase/Diguanylate cyclase domain", "RT-DGC", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Mononegavirales mRNA-capping domain V", "mRNA-CDV", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Hepatitis E virus, hinge domain", "X-domain", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Hepatitis E virus structural protein 2", "HEV-SP2", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("P/V phosphoprotein, paramyxoviral", "P/V-PP", interpro_description))


selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("3a-like VP, coronavirus", "3a-like VP", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("NendoU domain, nidovirus", "NendoU", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Hepatitis E virus Orf2, capsid", "HEV-Orf2-Capsid", interpro_description)) %>%
  mutate(interpro_description = gsub("Coronavirus Orf3a/b", "CoV-Orf3a/b", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Zinc finger, CCCH-type", "CCCH-ZF", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Assembly protein", "Assem-Protein", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Capsid protein C, pestivirus", "C-Protein", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Putative virion membrane protein", "VMP", interpro_description))

selected_rows <- selected_rows %>%
  mutate(interpro_description = gsub("Sialidase superfamily", "Sial", interpro_description))

genes <- as.data.frame(unique(selected_rows$interpro_description))
