metabo_agePstatus = ["1-arachidonoylglycerophosphocholine" ,
                     "1-eicosatrienoylglycerophosphocholine",
                     "1-linoleoylglycerophosphocholine",
                     "21-hydroxypregnenolone disulfate",    
                     "alpha-hydroxyisovalerate",    
                     "gamma-glutamylglutamine",
                     "glutamine",
                     "glycerophosphorylcholine GPC",
                     "lysine",        
                     "pregnen-diol disulfate"]

metabo_treat = ["2-hydroxybutyrate AHB",
                "arachidonate 204n6",
                "eicosapentaenoate EPA 205n3",
                "stearidonate 184n3",
                "alpha-tocopherol",
                "andro steroid monosulfate 2",
                "dihomo-linolenate 203n3 or n6",
                "docosahexaenoate DHA 226n3",
                "eicosapentaenoate EPA 205n3",   
                "glycerol 3-phosphate G3P",
                "linolenate alpha or gamma 183n3 or 6"]

metabo_timePstatus = ["1-arachidonylglycerol",
                      "1-docosahexaenoylglycerol 1-monodocosahexaenoin",
                      "5-oxoproline",
                      "cortisol",
                      "eicosapentaenoate EPA 205n3",
                      "glycerol",
                      "linolenate alpha or gamma 183n3 or 6",
                      "margarate 170",
                      "myristate 140",
                      "myristoleate 141n5",
                      "nonadecanoate 190",
                      "pregnen-diol disulfate",
                      "1-palmitoylglycerophosphate"]


metabo_sex = ['palmitoylcarnitine']


shapes_dict={'triangle':metabo_agePstatus,
             'square':  metabo_treat,
             'diamond': metabo_timePstatus,
             '*': metabo_sex}

young_age_thresh = 61

green_mets = ['eicosapentaenoate EPA 205n3',
 'glycerol',
 'linolenate alpha or gamma 183n3 or 6',
 'margarate 170',
 'myristate 140',
 'myristoleate 141n5',
 'nonadecanoate 190',
 'arachidonate 204n6',
 'dihomo-linolenate 203n3 or n6',
 'docosahexaenoate DHA 226n3',
 'stearidonate 184n3',
 '1-arachidonoylglycerophosphocholine',
 '1-eicosatrienoylglycerophosphocholine',
 '1-linoleoylglycerophosphocholine',
 '1-arachidonoylglycerophosphoinositol',
 '1-docosahexaenoylglycerophosphocholine',
 '1-eicosadienoylglycerophosphocholine',
 '1-heptadecanoylglycerophosphocholine',
 '1-myristoylglycerophosphocholine',
 '1-oleoylglycerophosphocholine',
 '1-palmitoleoylglycerophosphocholine',
 '1-palmitoylglycerophosphocholine',
 '1-stearoylglycerophosphocholine',
 '10-heptadecenoate 171n7',
 '10-nonadecenoate 191n9',
 '13-HODE  9-HODE',
 '16-hydroxypalmitate',
 '17-methylstearate',
 '2-palmitoylglycerophosphocholine',
 '2-stearoylglycerophosphocholine',
 '5-dodecenoate 121n7',
 'adrenate 224n6',
 'arachidate 200',
 'caprate 100',
 'cis-vaccenate 181n7',
 'dihomo-linoleate 202n6',
 'docosapentaenoate n3 DPA 225n3',
 'docosapentaenoate n6 DPA 225n6',
 'eicosenoate 201n9 or 11',
 'laurate 120',
 'linoleate 182n6',
 'oleate 181n9',
 'palmitate 160',
 'palmitoleate 161n7',
 'pentadecanoate 150',
 'stearate 180']

group2 = ['3-methyl-2-oxovalerate',
 '4-methyl-2-oxopentanoate',
 'fructose',
 'alpha-ketobutyrate',
 'gluconate',
 'glucose',
 'mannose',
 'cortisone']

group3 = ['1-palmitoylglycerophosphoethanolamine',
 '1-arachidonoylglycerophosphoethanolamine',
 '1-linoleoylglycerophosphoethanolamine',
 '2-linoleoylglycerophosphoethanolamine',
 '2-palmitoylglycerophosphoethanolamine',
 '1-oleoylglycerophosphoethanolamine',
 '2-oleoylglycerophosphoethanolamine',
 '1-stearoylglycerophosphoethanolamine']

group4 = ['lysine',
 'phenylalanine',
 'gamma-glutamylmethionine',
 'tryptophan',
 'methionine',
 'tyrosine',
 'isoleucine',
 'valine',
 'leucine',
 'gamma-glutamylleucine',
 'gamma-glutamylvaline',
 'gamma-glutamylglutamate',
 'gamma-glutamylphenylalanine',
 'gamma-glutamyltyrosine',
 'arginine']

group5 = ['2-hydroxybutyrate AHB',
 'acetylcarnitine',
 '3-hydroxybutyrate BHBA',
 '3-methyl-2-oxobutyrate',
 'hexanoylcarnitine',
 'octanoylcarnitine',
 'decanoylcarnitine',
 'cis-4-decenoyl carnitine',
 '2-aminobutyrate',
 'bilirubin Z,Z',
 'alpha-ketobutyrate',
 'mannose',
 '10-undecenoate 111n1',
 'riboflavin Vitamin B2',
 '3-methyl-2-oxovalerate',
 '4-methyl-2-oxopentanoate']

omega_3_mets = ['eicosapentaenoate_EPA_205n3',
'linolenate_alpha_or_gamma_183n3_or_6',
'arachidonate_204n6',
'dihomo-linolenate_203n3_or_n6',
'docosahexaenoate_DHA_226n3',
'stearidonate_184n3',
'dihomo-linoleate_202n6',
'linoleate_182n6',
'docosapentaenoate_n6_DPA_225n6',
'docosapentaenoate_n3_DPA_225n3',
'adrenate_224n6'
]

omega_3_mets = [c.replace("_", " ") for c in omega_3_mets]