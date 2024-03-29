rm(list=ls())


library(RTCGA)
library(dplyr)
library(RTCGAToolbox)
library(TCGAbiolinks)
library(data.table)
library(RDocumentation)
library(ggvenn)

setwd("F:/Ko/PD-1_responders/fgsea")


########################################################################################################################################################################################################################################################################################
#######################################################################    fgsea: REACTOME results frome M1T4 vs M3T4 comparison     ############################################################################################################################################
#######################################################################    T cell                                 ############################################################################################################################################
########################################################################################################################################################################################################################################################################################

M1T4_pos_reactome = c("REACTOME_SIGNALING_BY_INTERLEUKINS",	"REACTOME_SARS_COV_2_HOST_INTERACTIONS",	"REACTOME_PTEN_REGULATION",	"REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",	"REACTOME_INTERLEUKIN_1_SIGNALING",	"REACTOME_DNA_REPLICATION_PRE_INITIATION",	"REACTOME_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING",	"REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION",	"REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",	"REACTOME_NERVOUS_SYSTEM_DEVELOPMENT",	"REACTOME_REGULATION_OF_PTEN_STABILITY_AND_ACTIVITY",	"REACTOME_DEGRADATION_OF_GLI1_BY_THE_PROTEASOME",	"REACTOME_HEDGEHOG_ON_STATE",	"REACTOME_DOWNREGULATION_OF_ERBB2_SIGNALING",	"REACTOME_SELECTIVE_AUTOPHAGY",	"REACTOME_UCH_PROTEINASES",	"REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT",	"REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY",	"REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",	"REACTOME_TRAF6_MEDIATED_INDUCTION_OF_TAK1_COMPLEX_WITHIN_TLR4_COMPLEX",	"REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",	"REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT",	"REACTOME_RUNX1_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_DIFFERENTIATION_OF_HSCS",	"REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1",	"REACTOME_DOWNREGULATION_OF_ERBB4_SIGNALING",	"REACTOME_CRISTAE_FORMATION",	"REACTOME_REGULATION_OF_MRNA_STABILITY_BY_PROTEINS_THAT_BIND_AU_RICH_ELEMENTS",	"REACTOME_HEDGEHOG_LIGAND_BIOGENESIS",	"REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",	"REACTOME_AMYLOID_FIBER_FORMATION",	"REACTOME_MITOPHAGY",	"REACTOME_ALPHA_PROTEIN_KINASE_1_SIGNALING_PATHWAY",	"REACTOME_NOTCH2_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",	"REACTOME_TRANSLATION",	"REACTOME_SCAVENGING_BY_CLASS_F_RECEPTORS",	"REACTOME_RHO_GTPASES_ACTIVATE_RHOTEKIN_AND_RHOPHILINS",	"REACTOME_DEGRADATION_OF_AXIN",	"REACTOME_DEFECTIVE_CFTR_CAUSES_CYSTIC_FIBROSIS",	"REACTOME_DEGRADATION_OF_DVL",	"REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",	"REACTOME_STABILIZATION_OF_P53",	"REACTOME_MITOCHONDRIAL_TRANSLATION",	"REACTOME_REGULATION_OF_RAS_BY_GAPS",	"REACTOME_ABC_TRANSPORTER_DISORDERS",	"REACTOME_ASYMMETRIC_LOCALIZATION_OF_PCP_PROTEINS",	"REACTOME_RRNA_PROCESSING",	"REACTOME_PCP_CE_PATHWAY",	"REACTOME_SIGNALING_BY_ROBO_RECEPTORS",	"REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND",	"REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY",	"REACTOME_MAPK6_MAPK4_SIGNALING",	"REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",	"REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",	"REACTOME_NONSENSE_MEDIATED_DECAY_NMD",	"REACTOME_INFLUENZA_INFECTION",	"REACTOME_AUF1_HNRNP_D0_BINDS_AND_DESTABILIZES_MRNA",	"REACTOME_COMPLEX_I_BIOGENESIS",	"REACTOME_CELLULAR_RESPONSE_TO_STARVATION",	"REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",	"REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",	"REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",	"REACTOME_SELENOAMINO_ACID_METABOLISM",	"REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",	"REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",	"REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",	"REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",	"REACTOME_ATTENUATION_PHASE",	"REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",	"REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS",	"REACTOME_HSF1_ACTIVATION",	"REACTOME_HSF1_DEPENDENT_TRANSACTIVATION"
)
M3T4_pos_reactome = c("REACTOME_NETRIN_1_SIGNALING"
)

Nonreponder_pos_reactome =c("REACTOME_INTERLEUKIN_18_SIGNALING"
)
#RespondervsNonreponder_neg_reactome = c("REACTOME_INTERLEUKIN_18_SIGNALING")



data_1 = M1T4_pos_reactome 
data_2 = M3T4_pos_reactome
data_3 = Nonreponder_pos_reactome 
#data_4 = RespondervsNonreponder_neg_reactome


set1 = as.list(data_1)
set2 = as.list(data_2)
set3 = as.list(data_3)
#set4 = as.list(data_4)

venn.diagram(x=list(set1, set2, set3), 
             category.names = c("MaT4","McT4","PBMC-NR"), filename = "venndiagram_Tcell_REACTOME_only_(MaT4_vs._McT4)_Non-Responder.png", output=T, 
             imagetype = "png", height = 2000, width = 2000, resolution = 300, compression = "lzw",
             fill=c("skyblue", "yellowgreen","pink"), 
             cex = 2,  lwd=2, lty='blank', fontface="bold", fontfamily="sans", cat.cex=1, cat.fontface="bold", cat.default.pos="outer", cat.pos=c(0, 0, 200),cat.fontfamily="sans")




a <- list("M1T4_pos_reactome" = c("REACTOME_SIGNALING_BY_INTERLEUKINS",	"REACTOME_SARS_COV_2_HOST_INTERACTIONS",	"REACTOME_PTEN_REGULATION",	"REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING",	"REACTOME_INTERLEUKIN_1_SIGNALING",	"REACTOME_DNA_REPLICATION_PRE_INITIATION",	"REACTOME_BETA_CATENIN_INDEPENDENT_WNT_SIGNALING",	"REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION",	"REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",	"REACTOME_NERVOUS_SYSTEM_DEVELOPMENT",	"REACTOME_REGULATION_OF_PTEN_STABILITY_AND_ACTIVITY",	"REACTOME_DEGRADATION_OF_GLI1_BY_THE_PROTEASOME",	"REACTOME_HEDGEHOG_ON_STATE",	"REACTOME_DOWNREGULATION_OF_ERBB2_SIGNALING",	"REACTOME_SELECTIVE_AUTOPHAGY",	"REACTOME_UCH_PROTEINASES",	"REACTOME_ABC_FAMILY_PROTEINS_MEDIATED_TRANSPORT",	"REACTOME_SARS_COV_2_MODULATES_HOST_TRANSLATION_MACHINERY",	"REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",	"REACTOME_TRAF6_MEDIATED_INDUCTION_OF_TAK1_COMPLEX_WITHIN_TLR4_COMPLEX",	"REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",	"REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT",	"REACTOME_RUNX1_REGULATES_TRANSCRIPTION_OF_GENES_INVOLVED_IN_DIFFERENTIATION_OF_HSCS",	"REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1",	"REACTOME_DOWNREGULATION_OF_ERBB4_SIGNALING",	"REACTOME_CRISTAE_FORMATION",	"REACTOME_REGULATION_OF_MRNA_STABILITY_BY_PROTEINS_THAT_BIND_AU_RICH_ELEMENTS",	"REACTOME_HEDGEHOG_LIGAND_BIOGENESIS",	"REACTOME_CELLULAR_RESPONSE_TO_CHEMICAL_STRESS",	"REACTOME_AMYLOID_FIBER_FORMATION",	"REACTOME_MITOPHAGY",	"REACTOME_ALPHA_PROTEIN_KINASE_1_SIGNALING_PATHWAY",	"REACTOME_NOTCH2_ACTIVATION_AND_TRANSMISSION_OF_SIGNAL_TO_THE_NUCLEUS",	"REACTOME_TRANSLATION",	"REACTOME_SCAVENGING_BY_CLASS_F_RECEPTORS",	"REACTOME_RHO_GTPASES_ACTIVATE_RHOTEKIN_AND_RHOPHILINS",	"REACTOME_DEGRADATION_OF_AXIN",	"REACTOME_DEFECTIVE_CFTR_CAUSES_CYSTIC_FIBROSIS",	"REACTOME_DEGRADATION_OF_DVL",	"REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",	"REACTOME_STABILIZATION_OF_P53",	"REACTOME_MITOCHONDRIAL_TRANSLATION",	"REACTOME_REGULATION_OF_RAS_BY_GAPS",	"REACTOME_ABC_TRANSPORTER_DISORDERS",	"REACTOME_ASYMMETRIC_LOCALIZATION_OF_PCP_PROTEINS",	"REACTOME_RRNA_PROCESSING",	"REACTOME_PCP_CE_PATHWAY",	"REACTOME_SIGNALING_BY_ROBO_RECEPTORS",	"REACTOME_HSP90_CHAPERONE_CYCLE_FOR_STEROID_HORMONE_RECEPTORS_SHR_IN_THE_PRESENCE_OF_LIGAND",	"REACTOME_CHAPERONE_MEDIATED_AUTOPHAGY",	"REACTOME_MAPK6_MAPK4_SIGNALING",	"REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",	"REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING",	"REACTOME_NONSENSE_MEDIATED_DECAY_NMD",	"REACTOME_INFLUENZA_INFECTION",	"REACTOME_AUF1_HNRNP_D0_BINDS_AND_DESTABILIZES_MRNA",	"REACTOME_COMPLEX_I_BIOGENESIS",	"REACTOME_CELLULAR_RESPONSE_TO_STARVATION",	"REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",	"REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",	"REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",	"REACTOME_SELENOAMINO_ACID_METABOLISM",	"REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",	"REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",	"REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT",	"REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",	"REACTOME_ATTENUATION_PHASE",	"REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",	"REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS",	"REACTOME_HSF1_ACTIVATION",	"REACTOME_HSF1_DEPENDENT_TRANSACTIVATION"
),
"M3T4_pos_reactome" = c("REACTOME_NETRIN_1_SIGNALING"),
"Nonreponder_pos_reactome" = c("REACTOME_INTERLEUKIN_18_SIGNALING")
#,
#"RespondervsNonreponder_neg_reactome" = c("REACTOME_INTERLEUKIN_18_SIGNALING")
)


#data_1 = M1T4vs.all_pos_reactome 
#data_2 = M1T4vs.all_neg_reactome
#data_3 = M2T3_pos_reactome 
#data_4 = M2T3_neg_reactome


ggvenn(a, c("M1T4_pos_reactome", "Nonreponder_pos_reactome"))            # draw two-set venn
ggvenn(a, c("M3T4_pos_reactome", "Nonreponder_pos_reactome"))
#ggvenn(a, c("Set 1", "Set 2", "Set 3"))   # draw three-set venn
ggvenn(a)   # without set names, the first 4 elements in list will be chose to draw four-set venn








M1T4pos_Nonresponderpos_intersect= as.data.frame(Reduce(intersect, list(set1,set3)))
M1T4pos_Nonresponderpos_intersect =t(M1T4pos_Nonresponderpos_intersect)
M3T4pos_Nonresponderpos_intersect= as.data.frame(Reduce(intersect, list(set2,set3)))
M3T4pos_Nonresponderpos_intersect =t(M3T4pos_Nonresponderpos_intersect)
M1T4_M3T4_Nonresponder_all_intsect = as.data.frame(Reduce(intersect, list(set1,set2,set3)))
M1T4_M3T4_Nonresponder_all_intsect =t(M1T4_M3T4_Nonresponder_all_intsect)
M3T4_Nonresponder_only <- as.data.frame(Reduce(setdiff,list(M3T4pos_Nonresponderpos_intersect, M1T4_M3T4_Nonresponder_all_intsect)))
M3T4_Nonresponder_only =t(M3T4_Nonresponder_only)
M1T4_Nonresponder_only <- as.data.frame(Reduce(setdiff,list(M1T4pos_Nonresponderpos_intersect, M1T4_M3T4_Nonresponder_all_intsect)))
M1T4_Nonresponder_only =t(M1T4_Nonresponder_only)



n<-max(length(M1T4pos_Nonresponderpos_intersect),length(M3T4pos_Nonresponderpos_intersect),length(M1T4_M3T4_Nonresponder_all_intsect),length(M3T4_Nonresponder_only),length(M1T4_Nonresponder_only))
length(M1T4pos_Nonresponderpos_intersect)<-n
length(M3T4pos_Nonresponderpos_intersect)<-n
length(M1T4_M3T4_Nonresponder_all_intsect)<-n
length(M3T4_Nonresponder_only)<-n
length(M1T4_Nonresponder_only)<-n


M1T4_M3T4_responder_intersect = as.data.frame(cbind(M1T4pos_Nonresponderpos_intersect,M3T4pos_Nonresponderpos_intersect,M1T4_M3T4_Nonresponder_all_intsect,M3T4_Nonresponder_only,M1T4_Nonresponder_only))
colnames(M1T4_M3T4_responder_intersect)

M1T4_M3T4_responder_intersect<- M1T4_M3T4_responder_intersect %>% rename("V1" = "M1T4pos_responderpos_intersect")


write.csv(M1T4_M3T4_responder_intersect, "_Tcell_only_(MaT4_vs_M3T4)_Non-responder_REACTOME_intersect.csv")


########################################################################################################################################################################################################################################################################################
#######################################################################    fgsea: GOBP results from M1T4 vs M3T4 comparison     ############################################################################################################################################
#######################################################################    T cell                                 ############################################################################################################################################
########################################################################################################################################################################################################################################################################################

M1T4_pos_gobp = c("GOBP_ORGANOPHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_MITOCHONDRION_ORGANIZATION",	"GOBP_REGULATION_OF_CYSTEINE_TYPE_ENDOPEPTIDASE_ACTIVITY",	"GOBP_REGULATION_OF_PEPTIDASE_ACTIVITY",	"GOBP_NUCLEOBASE_CONTAINING_SMALL_MOLECULE_METABOLIC_PROCESS",	"GOBP_MITOCHONDRIAL_GENE_EXPRESSION",	"GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",	"GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS",	"GOBP_PROTEIN_LOCALIZATION_TO_MITOCHONDRION",	"GOBP_POSITIVE_REGULATION_OF_NF_KAPPAB_TRANSCRIPTION_FACTOR_ACTIVITY",	"GOBP_MITOCHONDRIAL_TRANSLATION",	"GOBP_MITOCHONDRIAL_TRANSPORT",	"GOBP_NUCLEOSIDE_PHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN",	"GOBP_CYTOPLASMIC_TRANSLATION",	"GOBP_HYPOTHALAMUS_CELL_DIFFERENTIATION",	"GOBP_HYPOTHALAMUS_GONADOTROPHIN_RELEASING_HORMONE_NEURON_DIFFERENTIATION",	"GOBP_POSITIVE_REGULATION_OF_NUCLEOTIDE_BINDING_OLIGOMERIZATION_DOMAIN_CONTAINING_2_SIGNALING_PATHWAY",	"GOBP_POSITIVE_REGULATION_OF_NUCLEOTIDE_BINDING_OLIGOMERIZATION_DOMAIN_CONTAINING_SIGNALING_PATHWAY",	"GOBP_ENERGY_DERIVATION_BY_OXIDATION_OF_ORGANIC_COMPOUNDS",	"GOBP_RIBOSE_PHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_ACTIVATION_OF_CYSTEINE_TYPE_ENDOPEPTIDASE_ACTIVITY_INVOLVED_IN_APOPTOTIC_PROCESS",	"GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_POTENTIAL",	"GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",	"GOBP_PURINE_CONTAINING_COMPOUND_BIOSYNTHETIC_PROCESS",	"GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",	"GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_1_PRODUCTION",	"GOBP_NADH_DEHYDROGENASE_COMPLEX_ASSEMBLY",	"GOBP_PROTON_TRANSMEMBRANE_TRANSPORT",	"GOBP_PROTEIN_LOCALIZATION_TO_PRESYNAPSE",	"GOBP_NUCLEOSIDE_TRIPHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS",	"GOBP_DEFENSE_RESPONSE_TO_GRAM_POSITIVE_BACTERIUM",	"GOBP_RETINA_HOMEOSTASIS",	"GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE",	"GOBP_ATP_METABOLIC_PROCESS",	"GOBP_ELECTRON_TRANSPORT_CHAIN",	"GOBP_CHAPERONE_MEDIATED_PROTEIN_COMPLEX_ASSEMBLY",	"GOBP_CELLULAR_RESPIRATION",	"GOBP_SKELETAL_MUSCLE_CONTRACTION",	"GOBP_DEFENSE_RESPONSE_TO_GRAM_NEGATIVE_BACTERIUM",	"GOBP_NEGATIVE_REGULATION_OF_INCLUSION_BODY_ASSEMBLY",	"GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE",	"GOBP_VASODILATION",	"GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",	"GOBP_MICROTUBULE_BASED_PROTEIN_TRANSPORT",	"GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_8_PRODUCTION",	"GOBP_ATP_BIOSYNTHETIC_PROCESS",	"GOBP_POSITIVE_REGULATION_OF_MUSCLE_CONTRACTION",	"GOBP_REGULATION_OF_INCLUSION_BODY_ASSEMBLY",	"GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",	"GOBP_RESPONSE_TO_TEMPERATURE_STIMULUS",	"GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN",	"GOBP_AEROBIC_RESPIRATION",	"GOBP_PROTEIN_FOLDING",	"GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS",	"GOBP_CHAPERONE_COFACTOR_DEPENDENT_PROTEIN_REFOLDING",	"GOBP_INCLUSION_BODY_ASSEMBLY",	"GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",	"GOBP_PROTEIN_REFOLDING",	"GOBP_CELLULAR_RESPONSE_TO_HEAT",	"GOBP_DE_NOVO_PROTEIN_FOLDING",	"GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING",	"GOBP_OXIDATIVE_PHOSPHORYLATION",	"GOBP_RESPONSE_TO_HEAT"
)
M3T4__pos_gobp = c("GOBP_REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS",	"GOBP_CAMP_METABOLIC_PROCESS"
)

Nonreponder_pos_gobp =c("GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL",	"GOBP_MULTICELLULAR_ORGANISMAL_IRON_ION_HOMEOSTASIS",	"GOBP_SYNAPTIC_MEMBRANE_ADHESION",	"GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL"
)
#RespondervsNonreponder_neg_reactome = c()



data_1 = M1T4_pos_gobp 
data_2 = M3T4__pos_gobp
data_3 = Nonreponder_pos_gobp 
#data_4 = RespondervsNonreponder_neg_reactome


set1 = as.list(data_1)
set2 = as.list(data_2)
set3 = as.list(data_3)
#set4 = as.list(data_4)

venn.diagram(x=list(set1, set2, set3), 
             category.names = c("MaT4","McT4","NR"), filename = "venndiagram_Tcell_GOBP_only_(MaT4_vs_McT4)_Nonresponder.png", output=T, 
             imagetype = "png", height = 2000, width = 2000, resolution = 300, compression = "lzw",
             fill=c("skyblue", "yellowgreen","pink"), 
             cex = 2,  lwd=2, lty='blank', fontface="bold", fontfamily="sans", cat.cex=1, cat.fontface="bold", cat.default.pos="outer", cat.pos=c(0, 0, 200),cat.fontfamily="sans")




a <- list("M1T4_pos_gobp" = c("GOBP_ORGANOPHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_MITOCHONDRION_ORGANIZATION",	"GOBP_REGULATION_OF_CYSTEINE_TYPE_ENDOPEPTIDASE_ACTIVITY",	"GOBP_REGULATION_OF_PEPTIDASE_ACTIVITY",	"GOBP_NUCLEOBASE_CONTAINING_SMALL_MOLECULE_METABOLIC_PROCESS",	"GOBP_MITOCHONDRIAL_GENE_EXPRESSION",	"GOBP_GENERATION_OF_PRECURSOR_METABOLITES_AND_ENERGY",	"GOBP_PURINE_CONTAINING_COMPOUND_METABOLIC_PROCESS",	"GOBP_PROTEIN_LOCALIZATION_TO_MITOCHONDRION",	"GOBP_POSITIVE_REGULATION_OF_NF_KAPPAB_TRANSCRIPTION_FACTOR_ACTIVITY",	"GOBP_MITOCHONDRIAL_TRANSLATION",	"GOBP_MITOCHONDRIAL_TRANSPORT",	"GOBP_NUCLEOSIDE_PHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN",	"GOBP_CYTOPLASMIC_TRANSLATION",	"GOBP_HYPOTHALAMUS_CELL_DIFFERENTIATION",	"GOBP_HYPOTHALAMUS_GONADOTROPHIN_RELEASING_HORMONE_NEURON_DIFFERENTIATION",	"GOBP_POSITIVE_REGULATION_OF_NUCLEOTIDE_BINDING_OLIGOMERIZATION_DOMAIN_CONTAINING_2_SIGNALING_PATHWAY",	"GOBP_POSITIVE_REGULATION_OF_NUCLEOTIDE_BINDING_OLIGOMERIZATION_DOMAIN_CONTAINING_SIGNALING_PATHWAY",	"GOBP_ENERGY_DERIVATION_BY_OXIDATION_OF_ORGANIC_COMPOUNDS",	"GOBP_RIBOSE_PHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_ACTIVATION_OF_CYSTEINE_TYPE_ENDOPEPTIDASE_ACTIVITY_INVOLVED_IN_APOPTOTIC_PROCESS",	"GOBP_REGULATION_OF_MITOCHONDRIAL_MEMBRANE_POTENTIAL",	"GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS",	"GOBP_PURINE_CONTAINING_COMPOUND_BIOSYNTHETIC_PROCESS",	"GOBP_MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",	"GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_1_PRODUCTION",	"GOBP_NADH_DEHYDROGENASE_COMPLEX_ASSEMBLY",	"GOBP_PROTON_TRANSMEMBRANE_TRANSPORT",	"GOBP_PROTEIN_LOCALIZATION_TO_PRESYNAPSE",	"GOBP_NUCLEOSIDE_TRIPHOSPHATE_BIOSYNTHETIC_PROCESS",	"GOBP_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS",	"GOBP_DEFENSE_RESPONSE_TO_GRAM_POSITIVE_BACTERIUM",	"GOBP_RETINA_HOMEOSTASIS",	"GOBP_POSITIVE_REGULATION_OF_NEUROINFLAMMATORY_RESPONSE",	"GOBP_ATP_METABOLIC_PROCESS",	"GOBP_ELECTRON_TRANSPORT_CHAIN",	"GOBP_CHAPERONE_MEDIATED_PROTEIN_COMPLEX_ASSEMBLY",	"GOBP_CELLULAR_RESPIRATION",	"GOBP_SKELETAL_MUSCLE_CONTRACTION",	"GOBP_DEFENSE_RESPONSE_TO_GRAM_NEGATIVE_BACTERIUM",	"GOBP_NEGATIVE_REGULATION_OF_INCLUSION_BODY_ASSEMBLY",	"GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_NADH_TO_UBIQUINONE",	"GOBP_VASODILATION",	"GOBP_RESPONSE_TO_TOPOLOGICALLY_INCORRECT_PROTEIN",	"GOBP_MICROTUBULE_BASED_PROTEIN_TRANSPORT",	"GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_8_PRODUCTION",	"GOBP_ATP_BIOSYNTHETIC_PROCESS",	"GOBP_POSITIVE_REGULATION_OF_MUSCLE_CONTRACTION",	"GOBP_REGULATION_OF_INCLUSION_BODY_ASSEMBLY",	"GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",	"GOBP_RESPONSE_TO_TEMPERATURE_STIMULUS",	"GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN",	"GOBP_AEROBIC_RESPIRATION",	"GOBP_PROTEIN_FOLDING",	"GOBP_PROTON_MOTIVE_FORCE_DRIVEN_ATP_SYNTHESIS",	"GOBP_CHAPERONE_COFACTOR_DEPENDENT_PROTEIN_REFOLDING",	"GOBP_INCLUSION_BODY_ASSEMBLY",	"GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",	"GOBP_PROTEIN_REFOLDING",	"GOBP_CELLULAR_RESPONSE_TO_HEAT",	"GOBP_DE_NOVO_PROTEIN_FOLDING",	"GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING",	"GOBP_OXIDATIVE_PHOSPHORYLATION",	"GOBP_RESPONSE_TO_HEAT"
),
"M3T4__pos_gobp" = c("GOBP_REGULATION_OF_EPITHELIAL_CELL_APOPTOTIC_PROCESS",	"GOBP_CAMP_METABOLIC_PROCESS"
),
"Nonreponder_pos_gobp" = c("GOBP_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL",	"GOBP_MULTICELLULAR_ORGANISMAL_IRON_ION_HOMEOSTASIS",	"GOBP_SYNAPTIC_MEMBRANE_ADHESION",	"GOBP_NATURAL_KILLER_CELL_MEDIATED_IMMUNE_RESPONSE_TO_TUMOR_CELL"
)
#,
#"RespondervsNonreponder_neg_reactome" = c("REACTOME_INTERLEUKIN_18_SIGNALING")
)


#data_1 = M1T4vs.all_pos_reactome 
#data_2 = M1T4vs.all_neg_reactome
#data_3 = M2T3_pos_reactome 
#data_4 = M2T3_neg_reactome


ggvenn(a, c("M1T4_pos_gobp", "Nonreponder_pos_gobp"))            # draw two-set venn
ggvenn(a, c("M3T4__pos_gobp", "Nonreponder_pos_gobp"))
#ggvenn(a, c("Set 1", "Set 2", "Set 3"))   # draw three-set venn
ggvenn(a)   # without set names, the first 4 elements in list will be chose to draw four-set venn








M1T4pos_Nonresponderpos_intersect= as.data.frame(Reduce(intersect, list(set1,set3)))
M1T4pos_Nonresponderpos_intersect =t(M1T4pos_Nonresponderpos_intersect)
M3T4pos_Nonresponderpos_intersect= as.data.frame(Reduce(intersect, list(set2,set3)))
M3T4pos_Nonresponderpos_intersect =t(M3T4pos_Nonresponderpos_intersect)
M1T4_M3T4_Nonresponder_all_intsect = as.data.frame(Reduce(intersect, list(set1,set2,set3)))
M1T4_M3T4_Nonresponder_all_intsect =t(M1T4_M3T4_Nonresponder_all_intsect)
M3T4_Nonresponder_only <- as.data.frame(Reduce(setdiff,list(M3T4pos_Nonresponderpos_intersect, M1T4_M3T4_Nonresponder_all_intsect)))
M3T4_Nonresponder_only =t(M3T4_Nonresponder_only)
M1T4_Nonresponder_only <- as.data.frame(Reduce(setdiff,list(M1T4pos_Nonresponderpos_intersect, M1T4_M3T4_Nonresponder_all_intsect)))
M1T4_Nonresponder_only =t(M1T4_Nonresponder_only)



n<-max(length(M1T4pos_Nonresponderpos_intersect),length(M3T4pos_Nonresponderpos_intersect),length(M1T4_M3T4_Nonresponder_all_intsect),length(M3T4_Nonresponder_only),length(M1T4_Nonresponder_only))
length(M1T4pos_Nonresponderpos_intersect)<-n
length(M3T4pos_Nonresponderpos_intersect)<-n
length(M1T4_M3T4_Nonresponder_all_intsect)<-n
length(M3T4_Nonresponder_only)<-n
length(M1T4_Nonresponder_only)<-n


M1T4_M3T4_responder_intersect = as.data.frame(cbind(M1T4pos_Nonresponderpos_intersect,M3T4pos_Nonresponderpos_intersect,M1T4_M3T4_Nonresponder_all_intsect,M3T4_Nonresponder_only,M1T4_Nonresponder_only))
colnames(M1T4_M3T4_responder_intersect)

M1T4_M3T4_responder_intersect<- M1T4_M3T4_responder_intersect %>% rename("V1" = "M1T4pos_responderpos_intersect")


write.csv(M1T4_M3T4_responder_intersect, "_Tcell_only_(MaT4_vs_M3T4)_Nonresponder_pos_GOBP_intersect.csv")




########################################################################################################################################################################################################################################################################################
#######################################################################    fgsea: KEGG results frome M1T4 vs M3T4 comparison     ############################################################################################################################################
#######################################################################    T cell                                 ############################################################################################################################################
########################################################################################################################################################################################################################################################################################

M1T4_pos_KEGG = c("KEGG_VIBRIO_CHOLERAE_INFECTION",	"KEGG_SPLICEOSOME",	"KEGG_RNA_POLYMERASE",	"KEGG_DRUG_METABOLISM_CYTOCHROME_P450",	"KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",	"KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",	"KEGG_PROTEIN_EXPORT",	"KEGG_HUNTINGTONS_DISEASE",	"KEGG_ALZHEIMERS_DISEASE",	"KEGG_CARDIAC_MUSCLE_CONTRACTION",	"KEGG_PARKINSONS_DISEASE",	"KEGG_OXIDATIVE_PHOSPHORYLATION",	"KEGG_RIBOSOME"
)
M3T4_pos_KEGG = c("KEGG_INTESTINAL_IMMUNE_NETWORK_FOR_IGA_PRODUCTION",	"KEGG_TYPE_I_DIABETES_MELLITUS",	"KEGG_ALLOGRAFT_REJECTION",	"KEGG_GRAFT_VERSUS_HOST_DISEASE",	"KEGG_AUTOIMMUNE_THYROID_DISEASE",	"KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",	"KEGG_PRIMARY_IMMUNODEFICIENCY"
)

Nonreponder_pos_KEGG =c(
)



data_1 = M1T4_pos_KEGG 
data_2 = M3T4_pos_KEGG
data_3 = Nonreponder_pos_KEGG 


set1 = as.list(data_1)
set2 = as.list(data_2)
set3 = as.list(data_3)

venn.diagram(x=list(set1, set2, set3), 
             category.names = c("MaT4","McT4","PBMC-NR"), filename = "venndiagram_Tcell_KEGG_only_(MaT4_vs._McT4)_Non-Responder.png", output=T, 
             imagetype = "png", height = 2000, width = 2000, resolution = 300, compression = "lzw",
             fill=c("skyblue", "yellowgreen","pink"), 
             cex = 2,  lwd=2, lty='blank', fontface="bold", fontfamily="sans", cat.cex=1, cat.fontface="bold", cat.default.pos="outer", cat.pos=c(0, 0, 200),cat.fontfamily="sans")




a <- list("M1T4_pos_KEGG" = c(),
"M3T4_pos_KEGG" = c(),
"Nonreponder_pos_KEGG" = c()
)


#data_1 = M1T4vs.all_pos_reactome 
#data_2 = M1T4vs.all_neg_reactome
#data_3 = M2T3_pos_reactome 
#data_4 = M2T3_neg_reactome


ggvenn(a, c("M1T4_pos_KEGG", "Nonreponder_pos_KEGG"))            # draw two-set venn
ggvenn(a, c("M3T4_pos_KEGG", "Nonreponder_pos_KEGG"))
#ggvenn(a, c("Set 1", "Set 2", "Set 3"))   # draw three-set venn
ggvenn(a)   # without set names, the first 4 elements in list will be chose to draw four-set venn








M1T4pos_Nonresponderpos_intersect= as.data.frame(Reduce(intersect, list(set1,set3)))
M1T4pos_Nonresponderpos_intersect =t(M1T4pos_Nonresponderpos_intersect)
M3T4pos_Nonresponderpos_intersect= as.data.frame(Reduce(intersect, list(set2,set3)))
M3T4pos_Nonresponderpos_intersect =t(M3T4pos_Nonresponderpos_intersect)
M1T4_M3T4_Nonresponder_all_intsect = as.data.frame(Reduce(intersect, list(set1,set2,set3)))
M1T4_M3T4_Nonresponder_all_intsect =t(M1T4_M3T4_Nonresponder_all_intsect)
M3T4_Nonresponder_only <- as.data.frame(Reduce(setdiff,list(M3T4pos_Nonresponderpos_intersect, M1T4_M3T4_Nonresponder_all_intsect)))
M3T4_Nonresponder_only =t(M3T4_Nonresponder_only)
M1T4_Nonresponder_only <- as.data.frame(Reduce(setdiff,list(M1T4pos_Nonresponderpos_intersect, M1T4_M3T4_Nonresponder_all_intsect)))
M1T4_Nonresponder_only =t(M1T4_Nonresponder_only)



n<-max(length(M1T4pos_Nonresponderpos_intersect),length(M3T4pos_Nonresponderpos_intersect),length(M1T4_M3T4_Nonresponder_all_intsect),length(M3T4_Nonresponder_only),length(M1T4_Nonresponder_only))
length(M1T4pos_Nonresponderpos_intersect)<-n
length(M3T4pos_Nonresponderpos_intersect)<-n
length(M1T4_M3T4_Nonresponder_all_intsect)<-n
length(M3T4_Nonresponder_only)<-n
length(M1T4_Nonresponder_only)<-n


M1T4_M3T4_responder_intersect = as.data.frame(cbind(M1T4pos_Nonresponderpos_intersect,M3T4pos_Nonresponderpos_intersect,M1T4_M3T4_Nonresponder_all_intsect,M3T4_Nonresponder_only,M1T4_Nonresponder_only))
colnames(M1T4_M3T4_responder_intersect)

M1T4_M3T4_responder_intersect<- M1T4_M3T4_responder_intersect %>% rename("V1" = "M1T4pos_Nonresponderpos_intersect")


write.csv(M1T4_M3T4_responder_intersect, "_Tcell_only_(MaT4_vs_M3T4)_Non-responder_KEGG_intersect.csv")

