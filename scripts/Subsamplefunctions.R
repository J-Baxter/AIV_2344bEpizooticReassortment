# Data Import and Subsampling Functions
FormatTipData <- function(x) {
  isolate_id <- x[grep("EPI_ISL_*", x)]
  
  subtype <- x[grep("H[[:digit:]]{1,}N[[:digit:]]*", x)] %>%
    str_extract(., "H[[:digit:]]{1,}N[[:digit:]]*")
  
  date <- x[grep("^[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]", x)]
  
  isolate_name <- x[grep("/[^/]*/[^/]*/", x)]
  
  if (identical(isolate_id, character(0))) {
    isolate_id <- NA
  }
  
  if (identical(subtype, character(0))) {
    subtype <- NA
  }
  
  if (identical(date, character(0))) {
    date <- NA
  }
  
  if (identical(isolate_name, character(0))) {
    isolate_name <- NA
  }
  
  out <- cbind.data.frame(
    "isolate_id" = isolate_id,
    "subtype" = subtype,
    "date" = date,
    "isolate_name" = isolate_name
  ) %>%
    distinct()
  
  
  return(out)
}


ExtractMetadata <- function(tiplabels){
  out <- tiplabels %>%
    gsub("\\|$", "", .) %>%
    sapply(., str_split_1, pattern = "\\|") %>%
    lapply(., FormatTipData) %>%
    bind_rows() %>%
    as_tibble() %>%
    
    # Format Dates
    mutate(date = parsedate::parse_date(date) %>%
             as.Date()) %>%
    mutate(decimal_date = format(round(decimal_date(date), 2), 
                                 nsmall = 2) ) %>%
    mutate(decimal_date = suppressWarnings(as.double(decimal_date))) %>%
    mutate(week_date = format(date, "%Y-%V")) %>%
    mutate(month_date = format(date, "%Y-%m")) %>%
    mutate(year_date = format(date, "%Y")) %>%
    
    # Format isolate name
    mutate(isolate_name = gsub("_(\\d{4})|(\\d{4})_", "\\1", perl = TRUE, isolate_name)) %>%
    mutate(isolate_name = case_when(!grepl("/(\\d{4})$", isolate_name) & !grepl("/$", isolate_name) ~ 
                                      paste0(isolate_name, "/NA"),
                                    !grepl("/(\\d{4})$", isolate_name) & grepl("/$", isolate_name) ~ paste0(isolate_name, "NA"),
                                    .default = isolate_name
    )) %>%
    
    # Extract data from isolate name
    separate_wider_delim(isolate_name,
                         delim = "/",
                         names = c("virus_species", "source", "location", "id_unsure", "year"),
                         cols_remove = F,
                         too_few = "align_end",
                         too_many = "debug"
    ) %>%
    
    # Drop extraneous cols and format
    dplyr::select(-c(year, contains("isolate_name_"))) %>%
    mutate(virus_species = case_when(str_length(source) == 1 ~ source,
                                     .default = virus_species
    )) %>%
    mutate(source = case_when(str_length(source) == 1 ~ NA,
                              .default = source
    )) %>%
    
    mutate(source = case_when(grepl("black[-_[:space:]]headed[-_[:space:]]heron|hartlaubs[-_[:space:]]gull|chicken|great[-_[:space:]]black[-_[:space:]]backed[-_[:space:]]gull|gallus[-_[:space:]]gallus|guinea[-_[:space:]]fowl|northern[-_[:space:]]shoveler", location, ignore.case = T) ~ location,
                              .default = source)) %>%
    mutate(location = case_when(grepl("black[-_[:space:]]headed[-_[:space:]]heron|hartlaubs[-_[:space:]]gull|chicken|great[-_[:space:]]black[-_[:space:]]backed[-_[:space:]]gull|gallus[-_[:space:]]gallus|guinea[-_[:space:]]fowl|northern[-_[:space:]]shoveler", location, ignore.case = T) & 
                                  grepl('africa|europe|kansas' , id_unsure) ~ id_unsure,
                                grepl("black[-_[:space:]]headed[-_[:space:]]heron|hartlaubs[-_[:space:]]gull|chicken|great[-_[:space:]]black[-_[:space:]]backed[-_[:space:]]gull|gallus[-_[:space:]]gallus|guinea[-_[:space:]]fowl|northern[-_[:space:]]shoveler", location, ignore.case = T) & 
                                  !grepl('africa|europe' , id_unsure) ~ NA,
                              .default = location))
  
  return(out)
}


FormatMetadata <- function(data){
  
  if('location_2' %in% colnames(data)){
    data <- data %>%
      # Format source before taxa allocation
      rename(joint_location = location) %>%
      mutate(across(contains('location'), 
                    .fns = ~ gsub("^NA$|^missing$|^unknown$", NA, .x))) %>%
      mutate(location = coalesce(location_3,location_2)) %>%
      select(-joint_location)%>%
      rename(source=host) %>%
      select(-starts_with('host')) %>%
      rename(month_date = date_year_month) %>%
      rename(decimal_date = date_frac) %>%
      mutate(year_date = date_year)
    
  } else{
    data <- data %>%
      mutate(date = parsedate::parse_date(date) %>%
               as.Date()) %>%
      mutate(decimal_date = format(round(decimal_date(date), 2), 
                                   nsmall = 2) ) %>%
      mutate(decimal_date = suppressWarnings(as.double(decimal_date))) %>%
      mutate(month_date = format(date, "%Y-%m")) %>%
      mutate(year_date = format(date, "%Y"))
  }
  
  data <- data %>%
    select(-matches('^date.+')) %>%
    #rename_with(.fn = ~gsub('_', '.', .x))%>%
    # Location (iso name and subdivision)
    mutate(location = tolower(location)) %>%
    mutate(location = gsub("[^A-Za-z]", " ", tolower(location))) %>%
    mutate(location =  gsub("\\s+", " ", str_trim(location))) %>%
    mutate(location = case_when(
      # Austria
      grepl('burgenland|stegersbach|mattersburg|bad sauerbrunn',
            location) ~ 'austria_burgenland',
      grepl('karnten', 
            location)  ~ 'austria_kärnten',
      grepl('niederosterreich|voesendorf|strasshof|stockerau|sollenau|breitensee|biberbach',
            location) ~ "austria_niederösterreich", 
      grepl('sankt poelten|pottschach|neunkirchen|markgrafneusiedl|koenigstetten|klosterneuburg',
            location) ~ "austria_niederösterreich",
      grepl('hof am leithaberge|haringsee|fishamend|bruderndorf|tulln|lower austria',
            location) ~ "austria_niederösterreich",
      grepl('oberosterreich$|sankt georgen|^linz$|oberneukirchen',
            location) ~ "austria_oberösterreich",
      grepl('salzburg', 
            location)  ~ 'austria_salzburg',
      grepl('steiermark|tagnitz|^graz$|eggersdorf|leibnitz', 
            location) ~ "austria_steiermark",
      grepl('tirol$', 
            location)  ~ 'austria_tirol',
      grepl('vorarlberg|voralrberg', 
            location)  ~ "austria_vorarlberg",
      grepl('^wien$|^vienna$|bundesland wien', 
            location) ~ 'austria_wien',
      grepl('carinthia', 
            location) ~ "austria_kärnten",
      
      # Belgium
      grepl("brussels capital region|region de bruxelles capitale",
            location) ~ "belgium_bruxelles",
      grepl('oost vlaanderen|west vlaandere|limburg|vlaams brabant|antwerp*',
            location) ~ "belgium_vlaanderen",
      grepl("brabant wallon|*hainaut$|*liege$| luxembourg|^namur$|province de namur",
            location) ~ "belgium_wallonie",
      
      
      # Bulgaria 
      grepl('blagoevgrad', 
            location) ~ "bulgaria_blagoevgrad",
      grepl('burgas$', 
            location) ~ "bulgaria_burgas",
      grepl('dobrich$', 
            location) ~ "bulgaria_dobrich",
      grepl('gabrovo$', 
            location) ~ "bulgaria_gabrovo",
      grepl('haskovo$', 
            location) ~ "bulgaria_haskovo",
      grepl('kardzhali$', 
            location) ~ "bulgaria_kardzhali",
      grepl('kyustendil$', 
            location) ~ "bulgaria_kyustendil",
      grepl('lovech$', 
            location) ~ "bulgaria_lovech",
      grepl('blagoevgrad', 
            location) ~ "bulgaria_blagoevgrad",
      grepl('"pazardzhik"', 
            location) ~ "bulgaria_pazardzhik",
      grepl('pernik', 
            location) ~ "bulgaria_pernik",
      grepl('pleven|slavyanovo', 
            location) ~ "bulgaria_pleven",
      grepl('plovdiv$', 
            location) ~ "bulgaria_plovdiv",
      grepl('^razgrad$', 
            location) ~ "bulgaria_razgrad",
      grepl('^ruse$', 
            location) ~ "bulgaria_ruse",
      grepl('shumen', 
            location) ~ "bulgaria_shumen",
      grepl('silistra', 
            location) ~ "bulgaria_silistra",
      grepl('sliven', 
            location) ~ "bulgaria_sliven",
      grepl('smolyan', 
            location) ~ "bulgaria_smolyan",
      grepl('^sofia$', 
            location) ~ "bulgaria_sofia",
      grepl('sofia-grad', 
            location) ~ "bulgaria_sofia-grad",
      grepl('stara zagora', 
            location) ~ "bulgaria_stara zagora",
      grepl('targovishte', 
            location) ~ "bulgaria_targovishte",
      grepl('varna', 
            location) ~ "bulgaria_varna",
      grepl('veliko tarnovo', 
            location) ~ "bulgaria_veliko tarnovo",
      grepl('^vidin$', 
            location) ~ "bulgaria_vidin",
      grepl('^vratsa$', 
            location) ~ "bulgaria_vratsa",
      grepl('^yambol$', 
            location) ~ "bulgaria_yambol",
      
      # Burkina Faso
      grepl("koubri", 
            location) ~ "burkina faso_kadiogo",
      grepl("gomboussougou", 
            location) ~ "burkina faso_zoundwéogo",
      
      # Cambodia
      grepl("krong phnum penh|phnom penh", 
            location) ~"cambodia_phnom penh",
      grepl("khett kandal",
            location) ~"cambodia_kândal",
      grepl("khett prey veng",
            location) ~"cambodia_prey vêng",
      grepl("takeo", 
            location) ~"cambodia_takêv",
      
      # Chile
      grepl("region de antofagasta", 
            location) ~ "chile_antofagasta",
      
      # Colombia
      grepl("departamento de bolivar", 
            location) ~"colombia_bolívar",
      grepl("departamento del choco", 
            location) ~"colombia_chocó",
      grepl("departamento de cordoba",
            location) ~"colombia_córdoba",
      grepl("departamento del magdalena", 
            location) ~"colombia_magdalena",
      
      # Costa Rica
      grepl("provincia de limon",
            location) ~ "costa rica_limón",
      
      # Croatia
      grepl('koprivnicko[ -]krizevacka', 
            location) ~ "croatia_koprivničko-križevačka",
      grepl('vukovarsko[- ]srijemska',
            location) ~ "croatia_vukovarsko-srijemska",
      grepl('sisacko[- ]moslavacka',
            location) ~ "croatia_sisacko-moslavacka",   
      
      # Denmark
      grepl('region sjaland', 
            location) ~    "denmark_sjælland",
      grepl('region syddanmark', 
            location) ~ "denmark_syddanmark",
      grepl('region nordjylland', 
            location) ~ "denmark_nordjylland",
      grepl('region midtjylland',
            location) ~"denmark_midtjylland",
      grepl('*denmark$',
            location) ~ 'denmark',
      
      # Ecuador
      grepl('provincia de manabi', 
            location) ~ "ecuador_manabi",
      
      # Honduras
      grepl('departamento de atlantida',
            location) ~ "honduras_atlántida",
      
      # Hungary
      grepl('bacs kiskun',
            location) ~   "hungary_bács-kiskun",
      grepl('csongrad megye', 
            location) ~ "hungary_csongrád",
      grepl('bekes megye', 
            location) ~ "hungary_békés",
      grepl('hajdu bihar', 
            location) ~ "hungary_hajdú-bihar",
      grepl('somogy megye', 
            location) ~  "hungary_somogy",
      
      # Iceland
      grepl('westfjords',
            location) ~ "iceland_vestfirðir",
      
      # Indonesia
      grepl('east java', 
            location) ~ "indonesia_jawa timur",
      
      # Iraq
      grepl('muhafazat ninawa',
            location) ~ "iraq_ninawa",
      
      # Israel
      grepl('northern district', 
            location) ~ "israel_hazafon",
      
      # Kazakhstan
      grepl('akmola', 
            location) ~ "kazakhstan_aqmola", 
      grepl('kostanay', 
            location) ~ "kazakhstan_qostanay",        
      grepl('soltüstik quzaqstan oblysy|north kazakhstan',
            location) ~ "kazakhstan_north kazakhstan",
      grepl('almaty province', 
            location) ~ "kazakhstan_almaty",
      
      # Latvia
      grepl('jurmala',
            location) ~"latvia_kurzeme",
      # Lithuania
      grepl('vilnius', 
            location) ~"lithuania_vilniaus",
      # Luxembourg
      grepl('district de grevenmacher', 
            location) ~ "luxembourg_grevenmacher",
      grepl('district de diekirch', 
            location) ~ "luxembourg_diekirch",
      # Mali
      grepl('kati', 
            location) ~ "mali_koulikoro",
      
      # Montenegro
      grepl('skadar lake',
            location) ~ "montenegro_bar",
      
      # Niger
      grepl('bouza tahoua', 
            location) ~ "niger_tahoua",
      grepl('torodi tillab ri|tillabéri',
            location) ~ "niger_tillabéry",  
      grepl(' zinder$', 
            location) ~ "niger_zinder",
      grepl('maradi$',
            location) ~ "niger_maradi",
      grepl('niamey$', 
            location) ~ "niger_niamey",
      
      # Norway
      grepl('hitra', 
            location) ~ "norway_sør-trøndelag",
      grepl('tromso',
            location) ~"norway_troms",
      
      # Peru
      grepl('departamento de lima', 
            location) ~  "peru_lima province",      
      grepl('provincia constitucional del callao',
            location) ~"peru_callao",
      grepl('departamento de tacna',
            location) ~ "peru_tacna",        
      grepl('departamento de ica',
            location) ~  "peru_ica",             
      
      # Portugal
      grepl('distrito de setubal', 
            location) ~ "portugal_setúbal",
      
      # Romania
      grepl('judetul timis', 
            location) ~ "romania_timiș",     
      grepl('judetul ilfov', 
            location) ~ "romania_ilfov",      
      grepl('judetul constanta', 
            location) ~ "romania_constanța", 
      grepl('mures ungheni', 
            location) ~ "romania_mureș",     
      grepl('judetul giurgiu', 
            location) ~ "romania_giurgiu",   
      grepl('tulcea',
            location) ~ "romania_tulcea",    
      
      # Senegal
      grepl('region de thies',
            location) ~ 'senegal_thiès',
      grepl('djoudj bird sanctuary',
            location) ~ 'sengeal_saint-louis',
      
      # Slovakia
      grepl('dobrohost|trnavský', 
            location) ~ "slovakia_trnavský",
      grepl('kosice|košický',
            location) ~  "slovakia_košický" ,
      grepl('kalinkovo|petrzalka|bratislavský',
            location) ~"slovakia_bratislavský",
      
      # South Africa
      grepl('province of the western cape',
            location) ~ "south africa_western cape",
      grepl('province of kwazulu natal', 
            location) ~  "south africa_kwazulu-natal",
      grepl('province of north west',
            location) ~ "south africa_north west",
      grepl('gauteng|pretoria|witwatersrand',
            location) ~ "south africa_gauteng",
      grepl('mpumalanga|eastern transvaal',
            location) ~ "south africa_gauteng",
      grepl('*free state$|vrystaat',
            location) ~ "south africa_free state",
      
      
      # Switzerland
      grepl('*zurich$',
            location) ~ "switzerland_zürich",
      grepl('*ticino$', 
            location) ~ "switzerland_ticino",    
      grepl('kanton basel stadt|*basel$', 
            location) ~ "switzerland_basel-stadt",
      grepl('*luzern$', 
            location) ~ "switzerland_lucerne",
      
      # Venezuela
      grepl("estado anzoategui", 
            location) ~ "venezuela_anzoátegui",
      
      # latin
      grepl('mexico', location) ~ 'méxico',
      grepl('^choco$', location) ~ 'colombia_chocó',
      grepl('valparaiso', location) ~ 'chile_valparaíso',
      grepl('tarapaca', location) ~ 'chile_tarapacá',
      grepl('^nuble$', location) ~ 'chile_ñuble',
      grepl('araucania', location) ~ 'chile_araucanía',
      #grepl("cordoba" , location) ~ 'spain_córdoba',
      grepl('atacama', location) ~ 'chile_atacama',
      grepl('coquimbo', location) ~ 'chile_coquimbo',
      grepl('antofagasta', location) ~ 'chile_antofagasta',
      grepl('maule', location) ~ 'chile_maule',
      grepl('arica y parinacota', location) ~ 'chile_arica y parinacota',
      grepl('bolivar', location) ~ 'ecuador_bolivar',
      grepl('magdalena', location) ~ 'colombia_magdalena',
      grepl('jalisco', location) ~ 'méxico_jalisco',
      
      #Oz
      grepl('south australia', location) ~ 'australia_south australia',
      grepl('victoria', location) ~ 'australia_victoria',
      
      
      # Vietnam
      grepl("viet nam", 
            location) ~ "vietnam",
      grepl("quang ninh" ,
            location) ~  "vietnam_quảng ninh" ,
      grepl('tinh lang son',
            location) ~ "vietnam_lạng sơn",
      grepl('tinh thanh hoa',
            location) ~ "vietnam_thanh hóa",
      grepl('nghe an|nhge an', 
            location) ~  "vietnam_nghệ an", 
      grepl('tinh quang tri',
            location) ~ "vietnam_quảng trị",
      grepl('nam ha tinh', 
            location) ~ "vietnam_hà tĩnh",
      grepl('thu do ha noi',
            location) ~ "vietnam_hà nội",
      
      # grepl('laos', location) ~ "lao people's democratic republic",
      .default = location)) %>%
    
    mutate(location = case_when(
      # PRC
      grepl("northern china|peking|beijing", 
            location) ~ "china_beijing", # Interpretation from paper - check with Lu
      grepl("chongqin*", 
            location) ~ "china_chongqing",
      grepl("shanghai",
            location) ~ "china_shanghai",
      grepl('tianjing',
            location) ~ 'china_tianjin',
      grepl('anhui',
            location) ~ 'china_anhui',
      grepl('^fujain$|^fujian$', 
            location) ~ 'china_fujian',
      grepl('gansu',
            location) ~ 'china_gansu',
      grepl('qingyuan|foshan|guangdong',
            location) ~ 'china_guangdong',
      grepl('guizhou', 
            location) ~ 'china_guizhou',
      grepl("heinan|hainan", 
            location) ~ "china_hainan", # confirm 
      grepl("wuhan|hebei", 
            location) ~ "china_hebei",
      grepl("heilongjiang", 
            location) ~ "china_heilongjiang",
      grepl("henan|sanmenxia|dongting", 
            location) ~ "china_henan",
      grepl("hubei",
            location) ~ "china_hubei",
      grepl("hunan|changsha", 
            location) ~ "china_hunan",
      grepl("jiangsu|suzhou|xuzhou|^xuyi$",
            location) ~ "china_jiangsu",
      grepl("jiangxi|ganzhou", 
            location) ~ "china_jiangxi",
      grepl("jilin", 
            location) ~ "china_jilin",
      grepl("liaoning",
            location) ~ "china_liaoning",
      grepl("qinghai", 
            location) ~ "china_qinghai",
      grepl("shaanxi", 
            location) ~ "china_shaanxi",
      grepl("shandong|rongcheng", 
            location) ~ "china_shandong",
      grepl("shanxi|southwestern china",
            location) ~ "china_shanxi",
      grepl("sichuan",
            location) ~ "china_sichuan",
      grepl("taiwan|chiayi|taoyuan|pingtung", 
            location) ~ "taiwan",
      grepl("taipei", 
            location) ~ "taiwan_taipei",
      grepl("yunnan.*", 
            location) ~ "china_yunnan",
      grepl("zhejiang|hangzhou|nanji$|eastern china",
            location) ~ "china_zhejiang",
      grepl("guangxi$|guilin|hechi|guangxi zhuang",
            location) ~ "china_guangxi",
      grepl("inner mongolia|tumuji|nei mongol", 
            location) ~ "china_nei mongol",
      grepl("ningxia$|ningxia hui", 
            location) ~ "china_ningxia hui",
      grepl("xinjiang$", 
            location) ~ "china_xinjiang uygur",
      grepl("xizang$|tibet", 
            location) ~ "china_xizang",
      grepl("macau",
            location) ~ "china_macau",
      grepl("hong kong|hongkong",
            location) ~ "china_hong kong",
      
      
      # Canada
      grepl("^ab$|alberta",
            location) ~ "canada_alberta",
      grepl("^bc$|british columbia", 
            location) ~ "canada_british columbia",
      grepl('^mb$|manitoba', 
            location) ~'canada_manitoba',
      grepl('^pei$|prince edward island',
            location) ~ 'canada_prince edward island',
      grepl('^qc$|quebec',
            location) ~ 'canada_québec',
      grepl('^on$|ontario',
            location) ~ 'canada_ontario',
      grepl('^sk$|saskatchewan',
            location) ~ 'canada_saskatchewan',
      grepl('new brunswick', 
            location) ~ 'canada_new brunswick',
      grepl('newfoundland*|*labrador$', 
            location) ~ 'canada_newfoundland and labrador',
      grepl('nova scotia$', 
            location) ~ 'canada_nova scotia',
      
      
      
      # Czechia
      grepl("jihocesky kraj|nov hrady by ov|české budějovice", 
            location) ~"czechia_jihočeský",
      grepl("south moravian region", 
            location) ~"czechia_jihomoravský",
      grepl("kralovehradecky kraj", 
            location) ~"czechia_královéhradecký",
      grepl("liberecky kraj",
            location) ~"czechia_liberecký",
      grepl("morovskoslezsky kraj", 
            location) ~"czechia_moravskoslezský",
      grepl("olomoucky kraj",
            location) ~"czechia_olomoucký",
      grepl("pardubicky kraj", 
            location) ~"czechia_pardubický",
      grepl("plzensky kraj",
            location) ~"czechia_plzeňský",
      grepl("hlavni mesto praha|prague", 
            location) ~"czechia_prague",
      grepl("stredocesky kraj",
            location) ~"czechia_středočeský",
      grepl("ustecky kraj",
            location) ~"czechia_ústecký",
      grepl("vysocina kraj|vysočina",
            location) ~"czechia_kraj vysočina",
      grepl("zlinsky kraj",
            location) ~"czechia_zlínský",
      
      
      # France
      grepl("^france&",
            location) ~"france",
      grepl("ain|allier|ardèche|cantal|drôme|isère|loire|haute-loire|puy-de-dôme|rhône|savoie|haute-savoie|rhone alpes", 
            location) ~ "france_auvergne-rhône-alpes",
      
      # bourgogne-franche-comté
      grepl("côte-d'or|doubs|jura|nièvre|haute-saône|saône-et-loire|yonne|territoire de belfort", 
            location) ~ "france_bourgogne-franche-comté",
      
      # brittany (bretagne)
      grepl("côtes-d'armor|finistère|ille-et-vilaine|morbihan|brittany",
            location) ~ "france_bretagne",
      
      # centre-val de loire
      grepl("^cher$|eure-et-loir|^indre$|indre-et-loire|loir-et-cher|loiret",
            location) ~ "centre-val de loire",
      
      # corsica (corse)
      grepl("corse-du-sud|haute-corse|^corsica$|^corse$", 
            location) ~ "france_corse",
      
      # grand est
      grepl("ardennes|^aube$|^marne$|haute-marne|meurthe-et-moselle|^meuse$|moselle|bas-rhin|haut-rhin|vosges|lorraine", 
            location) ~ "france_grand est",
      
      # hauts-de-france
      grepl("^aisne$|^nord$|^oise$|pas-de-calais|^somme$|nord pas de calais|picardie", 
            location) ~ "france_hauts-de-france",
      
      # île-de-france
      grepl("paris|seine-et-marne|yvelines|essonne|hauts-de-seine|seine-saint-denis|val-de-marne|val-d'oise|ile de france",
            location) ~ "france_île-de-france",
      
      # normandy (normandie)
      grepl("calvados|^eure$|manche|^orne$|seine-maritime|basse normandie", 
            location) ~ "france_normandie",
      
      # nouvelle-aquitaine
      grepl("charente|charente-maritime|corrèze|^creuse$|dordogne|gironde|landes|lot-et-garonne|pyrénées-atlantiques|deux-sèvres|vienne|haute-vienne|poitou charentes", 
            location) ~ "france_nouvelle-aquitaine",
      
      # occitanie
      grepl("ariège|^aude$|aveyron|^gard$|haute-garonne|^gers$|hérault|^lot$|lozère|hautes-pyrénées|pyrénées-orientales|tarn|tarn-et-garonne|midi pyrenees|languedoc roussillon", 
            location) ~ "france_occitanie",
      
      # pays de la loire
      grepl("loire-atlantique|maine-et-loire|mayenne|^sarthe$|vendée", 
            location) ~ "france_pays de la loire",
      
      # provence-alpes-côte d'azur
      grepl("alpes-de-haute-provence|hautes-alpes|alpes-maritimes|bouches-du-rhône|^var$|vaucluse", 
            location) ~ "france_provence-alpes-côte d'azur",
      
      
      # Italy
      grepl("l'aquila|teramo|pescara|^chieti$", 
            location) ~ "italy_abruzzo",
      grepl("potenza|matera", 
            location) ~ "italy_basilicata",
      grepl("catanzaro|cosenza|reggio calabria|crotone|vibo valentia", 
            location) ~ "italy_calabria",
      grepl("napoli|salerno|avellino|benevento|caserta", 
            location) ~ "italy_campania",
      
      grepl("bologna|ferrara|forlì-cesena|^modena$|^parma$|piacenza|ravenna|reggio emilia|^rimini$",
            location) ~ "italy_emilia-romagna",
      
      grepl("gorizia|pordenone|trieste|udine|^ud$", 
            location) ~ "italy_friuli-venezia giulia",
      
      grepl("roma|latina|rieti|viterbo|frosinone",
            location) ~ "italy_lazio",
      
      grepl("genova|imperia|la spezia|savona",
            location) ~ "italy_liguria",
      
      grepl("milano|bergamo|brescia|^bs$|^como$|^cr$|^cremona$|lecco|^lodi$|mantova|^mn$|pavia|sondrio|^varese$|lombardy", 
            location) ~ "italy_lombardia",
      
      grepl("ancona|ascoli piceno|fermo|macerata|pesaro e urbino", 
            location) ~ "italy_marche",
      
      grepl("campobasso|isernia", 
            location) ~ "italy_molise",
      
      grepl("torino|alessandria|^asti$|biella|^cuneo$|^novara$|verbano-cusio-ossola|vercelli", 
            location) ~ "italy_piemonte",
      
      grepl("bari|brindisi|foggia|lecce|taranto|barletta-andria-trani|puglia", 
            location) ~ "italy_apulia",
      
      grepl("cagliari|nuoro|oristano|sassari|sud sardegna", 
            location) ~ "italy_sardegna",
      
      grepl("agrigento|caltanissetta|catania|^enna$|messina|palermo|ragusa|siracusa|trapani",
            location) ~ "italy_sicily",
      
      grepl("firenze|arezzo|grosseto|livorno|lucca|massa-carrara|^pisa$|pistoia|prato|siena", 
            location) ~ "italy_toscana",
      
      grepl("trento|bolzano", 
            location) ~ "italy_trentino-alto adige",
      
      grepl("perugia|^terni$", 
            location) ~ "italy_umbria",
      
      grepl("^aosta$", 
            location) ~ "italy_valle d'aosta",
      
      grepl("belluno|padova|^pd$|rovigo|treviso|venezia|verona|vicenza|^v[eri]$|^veneto$", 
            location) ~ "italy_veneto",
      
      
      
      # Kosovo 
      grepl('kacanik|ferizaj|shtime|suhareke',
            location) ~ "kosovo_uroševac",
      grepl('gjilan',
            location) ~ "kosovo_gnjilane",  
      grepl('mitrovice|mitrovica', 
            location) ~ "kosovo_kosovska mitrovica",
      grepl('rahovec', 
            location) ~ "kosovo_đakovica",
      grepl('suhareke|malisheve',
            location) ~ "kosovo_prizren",
      grepl('kosovo', 
            location) ~ "kosovo",
      
      
      
      # Russian Federation
      grepl('russian federation', 
            location) ~ 'russia',
      grepl('khabarovsk', 
            location) ~ "russia_khabarovsk",
      grepl('chelyabinsk', 
            location) ~ "russia_chelyabinsk",
      grepl("astrakhan", 
            location) ~ "russia_astrakhan'",
      grepl("novosibirsk|chany lake|^chany$", 
            location) ~ "russia_novosibirsk",
      grepl("^omsk*|russia omsk region|russia omsk", 
            location) ~ "russia_omsk",
      grepl("tomsk*", 
            location) ~ "russia_tomsk",
      grepl("rostov oblast|rostov on don|^rostov$", 
            location) ~ "russia_rostov",
      grepl("russia primorje", 
            location) ~ "russia_primor'ye",
      grepl("sakhalin", 
            location) ~ "russia_sakhalin",
      grepl("yakutia", 
            location) ~ "russia_sakha",
      grepl("kostrom*", 
            location) ~ "russia_kostroma",
      grepl("amur region", 
            location) ~ "russia_amur",
      grepl("buryatia", 
            location) ~ "russia_buryat",
      grepl('north ossetia-alania', 
            location) ~ 'russia_north ossetia',
      grepl('dagestan', 
            location) ~ 'russia_dagestan',
      grepl('stavropol', 
            location) ~ "russia_stavropol'",
      grepl('krasnodar', 
            location) ~ "russia_krasnoyarsk", 
      grepl('tyumen|tumen', 
            location) ~ "russia_tyumen'",
      grepl("magadan",
            location) ~ "russia_magadan",
      grepl('saratov', 
            location) ~ "russia_saratov",
      grepl('kurgan*', 
            location) ~ "russia_kurgan",
      grepl('central russia',
            location) ~ 'russia',
      grepl('republic of tatarstan|^tatarstan$', 
            location) ~"russia_tatarstan",
      grepl('republic of kalmykia|^kalmykia$', 
            location) ~ "russia_kalmyk",
      grepl('siberia', 
            location) ~ "russia_krasnoyarsk",
      
      grepl('^tuva$|^tyva$|uvs nuur lake', 
            location) ~ "russia_tuva",
      
      grepl('kaliningra{0,1}d', 
            location) ~ "russia_kaliningrad",
      grepl('cheboksary', 
            location) ~ "russia_chuvash",
      grepl('kursk', 
            location) ~ "russia_kursk",
      grepl('sergiyev[[:space:]]{0,1}posad|shchyolkovo|moscow|chernogolovka', 
            location) ~ "russia_moskva",
      grepl('kamchatka', 
            location) ~ 'russia_kamchatka',
      grepl('voronezh', 
            location) ~ "russia_voronezh",
      grepl('^orel$', 
            location) ~ "russia_orel",
      grepl('^penza$', 
            location) ~ "russia_penza",
      grepl('^samara$', 
            location) ~ "russia_samara",
      grepl('^mari el$', 
            location) ~ 'russia_mariy-el',
      grepl('^leningrad$|leningrad region', 
            location) ~ "russia_leningrad",
      grepl('^orenburg$', 
            location) ~ "russia_orenburg",
      grepl('^ryazan$', 
            location) ~ "russia_ryazan'",
      grepl('^altai$', 
            location) ~ 'russia_altay',
      
      #morocco
      grepl('^fes$', 
            location) ~ 'morocco_fès - boulemane',
      grepl('^casablanca$', 
            location) ~ 'morocco_grand casablanca',
      grepl('^kenitra$', 
            location) ~ 'morocco_rabat - salé - zemmour - zaer',
      grepl('^nador$', 
            location) ~ 'morocco_oriental',
      
      
      # Egypt
      grepl('^luxor$', 
            location) ~ 'egypt_qina',
      grepl('^cairo$', 
            location) ~ 'egypt_al qahirah',
      grepl('al sharqia', 
            location) ~ 'egypt_ash sharqiyah',
      grepl('el wadialgidid', 
            location) ~ 'egypt_al wadi al jadid',
      
      
      # USA
      grepl("north dakota",
            location) ~ "united states_north dakota",
      grepl("massachussetts|conneticut|connecticut",
            location) ~"united states_connecticut",
      grepl("massachusetts",
            location) ~"united states_massachusetts",
      grepl("deleware|delaware", 
            location) ~"united states_delaware",
      grepl('^usa$', location) ~ 'united states',
      grepl('alaska$', location) ~ 'united states_alaska',
      grepl("lou[i]{0,1}siana", location) ~ 'united states_louisiana',
      grepl('alabama', location) ~ 'united states_alabama',
      grepl('arizona', location)~ 'united states_arizona',
      grepl('arkansas', location)~ 'united states_arkansas',
      grepl('california', location)~ 'united states_california',
      grepl('colorado', location)~ 'united states_colorado',
      grepl('district of columbia', location)~ 'united states_district of columbia',
      grepl('florida', location)~ 'united states_florida',
      grepl('georgia', location)~ 'united states_georgia',
      grepl('hawaii', location)~ 'united states_hawaii',
      grepl('idaho', location)~ 'united states_idaho',
      grepl('illinois', location)~ 'united states_illinois',
      grepl('indiana', location)~ 'united states_indiana',
      grepl('iowa', location)~ 'united states_iowa',
      grepl('kansas', location)~ 'united states_kansas',
      grepl('kentucky', location)~ 'united states_kentucky',
      grepl('maine', location)~ 'united states_maine',
      grepl( 'maryland', location)~ 'united states_maryland',
      grepl('massachusetts', location)~ 'united states_massachusetts',
      grepl('michigan', location)~ 'united states_michigan',
      grepl('minnesota', location)~ 'united states_minnesota',
      grepl('mississippi', location)~ 'united states_mississippi',
      grepl('missouri', location)~ 'united states_missouri',
      grepl('montana', location)~ 'united states_montana',
      grepl('nebraska', location)~ 'united states_nebraska',
      grepl('nevada', location)~ 'united states_nevada',
      grepl('new hampshire', location)~ 'united states_new hampshire',
      grepl('new jersey', location)~ 'united states_new jersey',
      grepl('new mexico', location)~ 'united states_new mexico',
      grepl('new york', location)~ 'united states_new york',
      grepl('north carolina', location)~ 'united states_north carolina',
      grepl('ohio', location)~ 'united states_ohio',
      grepl('oklahoma', location)~ 'united states_oklahoma',
      grepl('oregon', location)~ 'united states_oregon',
      grepl('pennsylvania', location)~ 'united states_pennsylvania',
      grepl('rhode island', location)~ 'united states_rhode island',
      grepl('south carolina', location)~ 'united states_south carolina',
      grepl('south dakota', location)~ 'united states_south dakota',
      grepl('tennessee', location)~ 'united states_tennessee',
      grepl('texas', location)~ 'united states_texas',
      grepl('utah', location)~ 'united states_utah',
      grepl('vermont', location)~ 'united states_vermont',
      grepl('virginia', location)~ 'united states_virginia',
      grepl('washington', location)~ 'united states_washington',
      grepl('west virginia', location)~ 'united states_west virginia',
      grepl('wisconsin', location)~ 'united states_wisconsin',
      grepl('wyoming', location)~ 'united states_wyoming',
      
      # UK
      grepl("wales|anglesey|county borough of wrexham",
            location) ~ "united kingdom_wales",
      grepl("england|barrow|northern ireland|county of kent|stratford|cheshire|derbyshire|northumberland|lincolnshire|dorset|yorkshire|hertfordshire|gloucestershire|*shire$|norfolk|essex|devon",
            location) ~ "united kingdom",
      grepl("scotland|orkney|western isles|fife|^angus$|dumfries*",
            location) ~ "united kingdom_scotland",
      
      # Japan
      grepl('aichi', 
            location) ~ 'japan_aichi',
      grepl('akita', 
            location) ~ 'japan_akita',
      grepl('aomori', 
            location) ~ 'japan_aomori',
      grepl('chiba', 
            location) ~ 'japan_chiba',
      grepl('ehime', 
            location) ~ 'japan_ehime',
      grepl('fukui', 
            location) ~ 'japan_fukui',
      grepl('fukuoka', 
            location) ~ 'japan_fukuoka',
      grepl('fukushima', 
            location) ~ 'japan_fukushima',
      grepl('gifu', 
            location) ~ 'japan_gifu',
      grepl('gunma', 
            location) ~ 'japan_gunma',
      grepl('hiroshima', 
            location) ~ 'japan_hiroshima',
      grepl('hokkaido', 
            location) ~ 'japan_hokkaido',
      grepl('hyogo|hyōgo', 
            location) ~ 'japan_hyōgo',
      grepl('tsukuba|ibaraki', 
            location) ~ 'japan_ibaraki',
      grepl('ishikawa', 
            location) ~ 'japan_ishikawa',
      grepl('iwate', 
            location) ~ 'japan_iwate',
      grepl('kagawa', 
            location) ~ 'japan_kagawa',
      grepl('kagoshima', 
            location) ~ 'japan_kagoshima',
      grepl('kanagawa', 
            location) ~ 'japan_kanagawa',
      grepl('kochi', 
            location) ~ 'japan_kochi',
      grepl('kumamoto', 
            location) ~ 'japan_kumamoto',
      grepl('kyoto', 
            location) ~ 'japan_kyoto',
      grepl('mie', 
            location) ~ 'japan_mie',
      grepl('miyagi', 
            location) ~ 'japan_miyagi',
      grepl('miyazaki', 
            location) ~ 'japan_miyazaki',
      grepl('nagano', 
            location) ~ 'japan_nagano',
      grepl('nagasaki|naoasaki', 
            location) ~ 'japan_naoasaki',
      grepl('nara', 
            location) ~ 'japan_nara',
      grepl('niigata', 
            location) ~ 'japan_niigata',
      grepl('oita', 
            location) ~ 'japan_oita',
      grepl('okayama', 
            location) ~ 'japan_okayama',
      grepl('okinawa', 
            location) ~ 'japan_okinawa',
      grepl('osaka', 
            location) ~ 'japan_osaka',
      grepl('saga', 
            location) ~ 'japan_saga',
      grepl('saitama', 
            location) ~ 'japan_saitama',
      grepl('shiga', 
            location) ~ 'japan_shiga',
      grepl('shimane', 
            location) ~ 'japan_shimane',
      grepl('shizuoka', 
            location) ~ 'japan_shizuoka',
      grepl('tochigi', 
            location) ~ 'japan_tochigi',
      grepl('tokushima', 
            location) ~ 'japan_tokushima',
      grepl('tokyo', 
            location) ~ 'japan_tokyo',
      grepl('tottori', 
            location) ~ 'japan_tottori',
      grepl('toyama', 
            location) ~ 'japan_toyama',
      grepl('wakayama', 
            location) ~ 'japan_wakayama',
      grepl('yamagata', 
            location) ~ 'japan_yamagata',
      grepl('yamaguchi', 
            location) ~ 'japan_yamaguchi',
      grepl('yamanashi', 
            location) ~ 'japan_yamanashi',
      # Spain CastillaLaMancha
      grepl("castilla[ ]{0,1}la[ ]{0,1}mancha|castille[ ]{0,1}la[ ]{0,1}mancha",
            location) ~ "spain_castilla-la mancha",
      
      # Poland
      grepl("lower silesian voivodeship|dolnośląskie", 
            location) ~"poland_dolnośląskie",
      grepl("kuyavian pomeranian voivodeship|kujawsko-pomorskie",
            location) ~"poland_kujawsko-pomorskie",
      grepl("lublin voivodeship", 
            location) ~"poland_lubelskie",
      grepl("lubusz voivodeship|lubuskie",
            location) ~"poland_lubuskie",
      grepl("lodzkie|odz voivodeship|łódzkie",
            location) ~"poland_łódzkie",
      grepl("lesser poland voivodeship|małopolskie", 
            location) ~"poland_małopolskie",
      grepl("masovian voivodeship|mazowieckie",
            location) ~"poland_mazowieckie",
      grepl("opole voivodeship|^opolskie$", 
            location) ~"poland_opolskie",
      grepl("pomeranian voivodeship|^pomorskie$",
            location) ~"poland_pomorskie",
      grepl("silesian voivodeship", 
            location) ~"poland_śląskie",
      grepl("swietokrzyskie voivodeship", 
            location) ~"poland_świętokrzyskie",
      grepl("warmi sko mazurskie",
            location) ~"poland_warmińsko-mazurskie",
      grepl("greater poland voivodeship|wielkopolskie", 
            location) ~"poland_wielkopolskie",
      grepl("west pomeranian voivodeship|zachodniopomorskie",
            location) ~"poland_zachodniopomorskie",
      
      
      # Sweden
      grepl("blekinge lan",
            location) ~"sweden_blekinge",
      grepl("dalarnas lan",
            location) ~"sweden_dalarna",
      grepl("gotlands lan",
            location) ~"sweden_gotland",
      grepl("hallands lan", 
            location) ~"sweden_halland",
      grepl("jonkopings lan", 
            location) ~"sweden_jönköping",
      grepl("kalmar lan",
            location) ~"sweden_kalmar",
      grepl("skane lan",
            location) ~"sweden_skåne",
      grepl("stockholms lan",
            location) ~"sweden_stockholm",
      grepl("sodermanlands lan", 
            location) ~"sweden_södermanland",
      grepl("uppsala lan", 
            location) ~"sweden_uppsala",
      grepl("vastra gotalands lan",
            location) ~"sweden_västra götaland",
      grepl("ostergotlands lan", 
            location) ~"sweden_östergötland",
      
      
      # Germany 
      grepl("^germany$", 
            location) ~"germany",
      grepl("germany[- ]bw|baden wuerttemberg|freiburg im breisgau|rosengarten|stutensee",
            location) ~"germany_baden-württemberg",
      grepl("germany[- ]by|bavaria|munich|nuremberg|wuerzburg|bayern", 
            location) ~"germany_bayern",
      grepl("germany[- ]hb|^bremen$", 
            location) ~"germany_bremen",
      grepl("germany[- ]hh|^hamburg$", 
            location) ~"germany_hamburg",
      grepl("germany[- ]he|hessen|freigericht|giessen|karben|ranstadt", 
            location) ~"germany_hessen",
      grepl("germany[- ]ni|lower saxony|brietlingen|dinklage|hannover|lueneburg|nordhorn|osnabruck",
            location) ~"germany_niedersachsen", 
      grepl("osterholz[- ]scharnbeck|theene|wingst",
            location) ~"germany_niedersachsen",
      grepl("germany[- ]nw|north rhine[- ]westphalia|aachen|brueggen", 
            location) ~"germany_nordrhein-westfalen",
      grepl("germany[- ]rp|rhineland[- ]palatinate", 
            location) ~"germany_rheinland-pfalz",
      grepl("germany[- ]sl|^saarland$", 
            location) ~"germany_saarland",
      grepl("germany[- ]sh|schleswig[- ]holstein|aumuehle|brickeln|luebeck|pinneberg|prohn", 
            location) ~"germany_schleswig-holstein",
      grepl("germany[- ]be|berlin",
            location) ~"germany_berlin",
      grepl("germany[- ]bb|brandenburg", 
            location) ~"germany_brandenburg",
      grepl("germany[- ]mv|mecklenburg[- ]vorpommern|grevesmuehlen|parchim|ruegen|warin", 
            location) ~"germany_mecklenburg-vorpommern",
      grepl("germany[- ]sn|sachsen|doberschuetz|dresden|gross dueben|leipzig|saxony",
            location) ~"germany_sachsen",
      grepl("germany[- ]st|sachsen[- ]anhalt|^halle$|zeitz", 
            location) ~"germany_sachsen-anhalt",
      grepl("germany[- ]th|thuring.*|^gera$", 
            location) ~"germany_thüringen",
      
      
      # Korea
      grepl("korea", 
            location) ~ "south korea",
      grepl('^gg$|gyeonggi*', 
            location) ~ "south korea_gyeonggi-do",
      grepl('^gw$', 
            location) ~ "south korea_gangwon-do",
      grepl('^cn$|chungcheongnam*', 
            location) ~ "south korea_chungcheongnam-do",
      grepl('^cb$|chungcheongbuk*',
            location) ~ "south korea_chungcheongbuk-do",
      grepl('^gn$',
            location) ~ "south korea_gyeongsangnam-do",
      grepl('^gb$',
            location) ~ "south korea_gyeongsangbuk-do",
      grepl('^jn$',
            location) ~ "south korea_jeollanam-do",
      grepl('^jb$', 
            location) ~ "south korea_jeollabuk-do",
      grepl('^jj$',
            location) ~ "south korea_jeju",
      grepl('^sw$', 
            location) ~ "south korea_seoul", 
      grepl('^ic$',
            location) ~ "south korea_incheon", 
      grepl('^dj$',
            location) ~ "south korea_daejeoni", 
      grepl('^sj$|sejong teugbyeolsi', 
            location) ~ "south korea_sejong", 
      grepl('^gj$',
            location) ~ "south korea_gwangju gwang", 
      grepl('^dg$',
            location) ~ "south korea_daegu", 
      grepl('^ps$', 
            location) ~ "south korea_busan", 
      grepl('^us$', 
            location) ~ "south korea_ulsan", 
      
      
      # Netherlands
      grepl('drenthe|westerveld|borger-odoorn|tynaarlo|eext|koekange|middendrexhe|noordenveld|nl assen',
            location)  ~ 'netherlands_drenthe',
      grepl('flevoland|almere|lelystad|noordoostpolder|noordosterpolder|marker wadden|stroe waddenkust',
            location) ~ 'netherlands_flevoland',
      grepl('zeewolde tulpeiland|zuidlaarmeergebiede oost polder|roggebotsluis|zeewolde|almeerder zand',
            location) ~ 'netherlands_flevoland',
      grepl('friesland|fryslan|tytsjerksteradiel|eastermar|heerenveen|schiermonnikoog|^ameland$', 
            location) ~ "netherlands_fryslân",
      grepl('kollum|^vlieland.*|ameland|noord[[:space:]]{0,1}buren waddenkust|greonterp|walterswald|nl abbega|terschelling|leeuwarden', 
            location) ~ "netherlands_fryslân",
      grepl('gelderland|zutphen|zevenaar|westvoort|westervoort|andelst|arnhem|bennekom|deliemers',
            location) ~ 'netherlands_gelderland',
      grepl('spijk|doetichem|doetinchem|epe|ermelo|gelmonde|ingen|klarenbeek|lochem|oldenbroek|huissen',
            location) ~ 'netherlands_gelderland',
      grepl('putten|zwolle bergkloosterweg|apeldoorn boogschutterstraat|boven leeuwen|barneveld',
            location) ~ 'netherlands_gelderland',
      grepl('groningen|zuidhorn|oldambt|spijk|lauwersmeer|haren boeremapark|scheemda|sint annen', 
            location)  ~ 'netherlands_groningen',
      grepl('limburg|venlo|gennep|kerkrade|landgraaf|ottersum|maastricht',
            location) ~ 'netherlands_limburg',
      grepl('noord brabant|wernhout|vlijmen|beekendonk|best|boekel|uden|agatha|schaik|rosmalen|^reek&',
            location) ~ 'netherlands_noord-brabant',
      grepl('lierop|oss|north brabant|oirschot|lithse eendenkooi|werkendam',
            location) ~ 'netherlands_noord-brabant',
      grepl('noord holland|bloemendaal|heerhugowaard|hilversum|huizen|naarden|reek|texel de schorren|volendam|durgerdam', 
            location) ~ 'netherlands_noord-holland',
      grepl('wie[r]{0,1}ingerwerf|den oever|hippolytushoef|kreupel|onderdijk|^obdam$|monnickendam|slootdorp|zuidoost beemster', 
            location) ~ 'netherlands_noord-holland',
      grepl('ijmuiden forteiland|zwaagdijk.*|schellinkhout zuiderdijk|schagen schagerweg|nl wormer$|west graftdijk|oostwoud', 
            location) ~ 'netherlands_noord-holland',
      grepl('overijssel|wierden|enschede|hardenberg|heino|losser|overdinkel|raalte|ijsselmuiden|nl zwolle|kamperveen|mastenbroek',
            location)  ~ 'netherlands_overijssel',
      grepl('utrecht|bilthoven|bosch en duin|bunnik|debilt|derondevenen|soest|haarzuilens|houten|rhenen', 
            location)  ~ 'netherlands_utrecht',
      grepl('ijsselstein|langbroek|nieuwegein|ijsselstein|woerden|leusden|vleuten.*|portengen|eemmeer dode hond|vianen', 
            location)  ~ 'netherlands_utrecht',
      grepl('zeeland|terneuzen|reimerswaal|grenspad|middelburg|weversinlaag|weversinlaag|neeltje jans', 
            location) ~ 'netherlands_zeeland',
      grepl('walcheren t vroon westkapelle|sloehaven', 
            location) ~ 'netherlands_zeeland',
      grepl('zuid holland|zoeterwoude|westland|den haag|rotterdam|rijswijk|leiden|oud albas|oud alblas',
            location) ~ 'netherlands_zuid-holland',
      grepl('nederlek|hardinxveldgiessendam|south holland|eendenkooi stolwijk|haringvliet|bliek$', 
            location) ~ 'netherlands_zuid-holland',
      grepl('kwade hoek stellendam|stellendam scheelhoekeiland|leidschendam|nl gouda|groot ammers|reeuwijk|stolwijk', 
            location) ~ 'netherlands_zuid-holland',
      grepl('^nl$|nl sovon', location) ~ 'netherlands',
      
      
      # Nigeria
      grepl("abuja",
            location) ~"nigeria_federal capital territory",
      grepl("abia state", 
            location) ~"nigeria_abia",
      grepl("adamawa state",
            location) ~"nigeria_adamawa",
      grepl("akwa ibom state", 
            location) ~"nigeria_akwa ibom",
      grepl("aguata lga|anambra state",
            location) ~"nigeria_anambra",
      grepl("tilden fulani toro lga|bauchi state", 
            location) ~"nigeria_bauchi",
      grepl("bayelsa state", 
            location) ~"nigeria_bayelsa",
      grepl("benue state", 
            location) ~"nigeria_benue",
      grepl("borno state", 
            location) ~"nigeria_borno",
      grepl("corss river state", 
            location) ~"nigeria_cross river",
      grepl("delta state",
            location) ~"nigeria_delta",
      grepl("ebonyi state", 
            location) ~"nigeria_ebonyi",
      grepl("benin edo|edo state", 
            location) ~"nigeria_edo",
      grepl("ekiti state",
            location) ~"nigeria_ekiti",
      grepl("enugu state", 
            location) ~"nigeria_enugu",
      grepl("gombe state|pantaini lbm|^akko[ gombe]{0,1}", 
            location) ~"nigeria_gombe",
      grepl("imo state", 
            location) ~"nigeria_imo",
      grepl("dutse lga|jigawa state",
            location) ~"nigeria_jigawa",
      grepl("goni gora chikun|kaduna state", 
            location) ~"nigeria_kaduna",
      grepl("kano state|gwale lga|gunduwawa gezawa|dawakin tofa|^kano$",
            location) ~"nigeria_kano",
      grepl("katsina state", 
            location) ~"nigeria_katsina",
      grepl("kebbi state",
            location) ~"nigeria_kebbi",
      grepl("kogi state", 
            location) ~"nigeria_kogi",
      grepl("kwara state",
            location) ~"nigeria_kwara",
      grepl("lagos state", 
            location) ~"nigeria_lagos",
      grepl("keffi lgc|nasarawa state",
            location) ~"nigeria_nasarawa",
      grepl("niger state", 
            location) ~"nigeria_niger",
      grepl("ogun state",
            location) ~"nigeria_ogun",
      grepl("ondo state", 
            location) ~"nigeria_ondo",
      grepl("osun state", 
            location) ~"nigeria_osun",
      grepl("oyo state", 
            location) ~"nigeria_oyo",
      grepl("plateau state|jos north|maraban jos|rayfield jos south|^plateau$", 
            location) ~"nigeria_plateau",
      grepl("port harcourt|rivers state|nkpor obio nkpor lgc", 
            location) ~"nigeria_rivers",
      grepl("sokoto state", 
            location) ~"nigeria_sokoto",
      grepl("taraba state",
            location) ~"nigeria_taraba",
      grepl("yobe state",
            location) ~"nigeria_yobe",
      grepl("zamfara state", 
            location) ~"nigeria_zamfara",
      
      grepl('bosnia - herzegovina',
            location) ~ 'bosnia and herzegovina',
      
      grepl('sva |^sva$|north america|zeebrugge belgi',
            location) ~ NA,
      
      # Indonesia
      grepl("hulu sungai utara|banjarbaru", 
            location) ~ "indonesia_kalimantan selatan",
      
      #SA
      grepl("al qasim", 
            location) ~ "saudi arabia_al qassim",
      grepl("riyadh", 
            location) ~ "saudi arabia_ar riyad",
      
      # Other
      grepl("republic of georgia", location) ~ "georgia",
      grepl("czech republic|czechia|czech reublic", location) ~ "czechia",
      grepl("burkina faso|burkina-faso", location) ~ "burkina faso",
      grepl('macedonia', location) ~ 'north macedonia',
      grepl('king island', location) ~ 'australia_tasmania',
      grepl('lisbon', location) ~ 'portugal_lisboa',
      grepl('ashanti', location) ~ 'ghana_ashanti',
      grepl('kerala', location) ~ 'india_kerala',
      grepl('^assam$', location) ~ 'india_assam',
      grepl('serres', location) ~ 'greece_macedonia and thrace',
      grepl('^south$', location) ~ 'iceland_suðurland',
      
      
      .default = location
    )) %>% 
    
    left_join(geodata, by = join_by(location == match)) %>%
    
    mutate(source = gsub("_", " ", tolower(source))) %>%
    mutate(source = str_trim(source)) %>%
    mutate(source= gsub('gray', 'grey', source)) %>%
    rowwise() %>%
    mutate(primary_com_name = FormatBird(source)) %>%
    mutate(primary_com_name = FormatMammal(primary_com_name)) %>%
    as_tibble() %>%
    mutate(primary_com_name = case_when(
      
      # location-specific birds
      grepl('sparrowhawk', 
            primary_com_name) & grepl('europe', region) ~ 'eurasian sparrowhawk',
      grepl("^magpie$|^common magpie",
            primary_com_name) & grepl('europe', region) ~ "eurasian magpie",
      grepl("^magpie$|^common magpie",
            primary_com_name) & grepl('eastern asia', region) ~ "oriental/eurasian magpie",
      grepl("^magpie$|^common magpie", 
            primary_com_name) & grepl('northern america', region) ~ "black-billed magpie",
      grepl("^kestrel$",
            primary_com_name) & grepl('northern america', region) ~ "american kestrel",
      grepl("^kestrel$", 
            primary_com_name) & grepl('europe', region) ~ "lesser/eurasian kestrel",
      grepl("^guinea {0,1}fowl", 
            primary_com_name) & grepl('africa', region) ~ "crested guineafowl sp.",
      grepl('^guinea {0,1}fowl', 
            primary_com_name) & !grepl('africa', region) ~ "helmeted guineafowl (domestic type)",
      grepl('^fulmar$', 
            primary_com_name) & !grepl('south america', region) ~ "northern fulmar",
      grepl('^fulmar$',
            primary_com_name) & grepl('south america', region) ~ "southern fulmar",
      grepl('^gannet$', 
            primary_com_name) & grepl('europe', region) ~ "northern gannet",
      grepl('jungle crow', 
            primary_com_name) ~ "crow/raven sp.",
      grepl('oystercatcher', 
            primary_com_name) & grepl('europe', region) ~ 'eurasian oystercatcher',
      grepl("turkey|^pavo$", 
            primary_com_name) & grepl('central america', region) ~ "wild turkey", 
      grepl("turkey|^pavo$",
            primary_com_name) & !grepl('central america', region) ~ "wild turkey (domestic type)", 
      grepl("turnstone", 
            primary_com_name) & grepl('europe', region) ~ "ruddy turnstone",
      primary_com_name == 'sea eagle' & grepl('europe', region) ~"white-tailed eagle",
      grepl("quail", 
            primary_com_name) & !grepl('america', region) ~ "old world quail sp.", 
      grepl("quail", 
            primary_com_name) & grepl('america', region) ~ "new world quail sp.", 
      grepl('^vulture$', 
            primary_com_name) &  grepl('america', region) ~ "new world vulture sp.",
      grepl('^vulture$', 
            primary_com_name) & !grepl('america', region) ~ "old world vulture sp.",
      
      # unknown
      grepl('^an$|unknown', 
            primary_com_name) ~ 'unknown',
      is.na(primary_com_name) ~ 'unknown',
      
      
      .default = primary_com_name
    )) %>%
    
    # Join taxonomy csvs
    left_join(birds, by = join_by(primary_com_name)) %>%
    left_join(mammals,  by = join_by(primary_com_name)) %>%
    
    # fix joins
    mutate(across(ends_with(".x"), .names = "{gsub('.x', '', {.col})}", 
                  .fns = ~coalesce(., get(sub("\\.x$", ".y", cur_column()))))) %>%
    select(-(ends_with('.x') |ends_with('.y'))) %>%
    
    # environment
    mutate(across(c(order, class, family),
                  ~ case_when(primary_com_name == 'environment' ~ 'environment', TRUE ~ .)),
           across(c(order, class, family),
                  ~ case_when(primary_com_name == 'unknown' ~ 'unknown', TRUE ~ .))) %>%
    
    # binary host information
    mutate(is_domestic = case_when(grepl("domestic", sci_name) ~ "domestic",
                                   grepl('homo sapiens',sci_name) ~ 'unknown',
                                   .default = "wild"
    )) %>%
    rowwise() %>%
    mutate(is_bird = isBird(tolower(order)))%>% 
    as_tibble() %>%
    
    # Format column names
    #dplyr::select(-c(virus_species, id_unsure)) %>%
    dplyr::rename(
      virus_subtype = subtype,
      collection_date = date,
      collection_datedecimal = decimal_date,
      collection_datemonth = month_date,
      collection_dateyear = year_date,
      collection_regionname = region,
      collection_countryname = country,
      collection_countrycode = gid_0,
      collection_countrylat = adm0_lat,
      collection_countrylong = adm0_long,
      collection_subdiv1name = name_1,
      collection_subdiv1code = hasc_1,
      collection_subdiv1lat = adm1_lat,
      collection_subdiv1long = adm1_long,
      #collection_subdiv2.name = name_subdiv2,
      #collection_subdiv2.code = code_subdiv2,
      #collection_subdiv2.lat = lat_subdiv2,
      #collection_subdiv2.long = long_subdiv2,
      host_class = class,
      host_order = order,
      host_family = family,
      host_sciname = sci_name,
      host_commonname = primary_com_name,
      host_isbird = is_bird,
      host_isdomestic = is_domestic
    ) %>%
    
    # Infer simplified host categories for BEAST
    mutate(host_simplifiedhost = case_when(
      host_order %in% c('anseriformes', 'galliformes', 'charadriiformes') ~ paste(host_order, 
                                                                                  host_isdomestic,
                                                                                  sep = '-'),
      host_class == 'mammalia' ~ 'mammal',
      host_commonname == 'unknown' ~ 'unknown',
      host_commonname == 'environment' ~ 'environment',
      .default = 'other')) %>%
    
    dplyr::relocate(
      virus_subtype,
      isolate_id,
      isolate_name,
      collection_date,
      collection_datedecimal,
      collection_datemonth,
      collection_dateyear,
      host_class,
      host_order,
      host_family,
      host_sciname,
      host_commonname,
      host_isbird,
      host_isdomestic,
      host_simplifiedhost,
      collection_regionname,
      collection_countryname,
      collection_countrycode,
      collection_countrylat,
      collection_countrylong,
      collection_subdiv1name,
      collection_subdiv1code,
      collection_subdiv1lat,
      collection_subdiv1long#,
      #collection_subdiv2_name,
      #collection_subdiv2_code,
      #collection_subdiv2_lat,
      #collection_subdiv2_long
    ) %>%
    
    # Select columns
    select(-c(source, location, varname_1, iso_1)) %>%
    
    # Format NAs
    mutate(across(everything(), .fns = ~ gsub('^NA$', NA, .x))) 
  
  
  return(data)
  
}


MergeReassortantData <-  function(data, newdata){
  df <- data  %>% 
    left_join(newdata, by = join_by(isolate_id)) %>%
    mutate(across(ends_with(".x"),
                  ~coalesce(., get(sub("\\.x$", ".y", cur_column()))),
                  .names = "{.col}")) %>%
    select(-ends_with(".y")) %>%
    rename_with(~gsub('.x', '', .x)) %>%
    select(-id_unsure) %>%
    mutate(collection_tipdate = case_when(
      is.na(collection_date) & is.na(collection_datemonth) ~ collection_dateyear,
      is.na(collection_date) & !is.na(collection_datemonth) ~ collection_datemonth,
      .default = as.character(collection_date))) 
  
  return(df)
}


# Where tiplabels are not in order, a neighbour joining tree formed and lables
ImputeCladeandCluster <- function(metadata, alignment, ordered = FALSE){
  
  seg <- metadata %>% 
    pull(segment) %>%
    unique()
  
  print(seg)
  
  if(ordered == FALSE){
    tree = nj(dist.dna(alignment))
    seqnames = tree$tip.label
    
  }else{
    seqnames = rownames(alignment)
  }
  
  out <- metadata %>%
    mutate(tipnames = ordered(tipnames,
                              levels = seqnames)) %>%
    arrange(tipnames) %>%
    
    # Impute clade
    mutate(clade = case_when(
      is.na(clade) & na.locf0(clade, fromLast = TRUE) == na.locf0(clade, fromLast = FALSE) ~ na.locf0(clade, fromLast = TRUE), 
      .default = clade)) %>%
    
    # Impute cluster (column is dependent on segment)
    rename(profile = cluster_profile) %>%
    rename(genome = cluster_genome) %>%
    pivot_longer(matches('cluster_[a-z]{2}\\d{0,1}$'), 
                 values_to = 'cluster_number',
                 names_to = 'cluster_segment') %>%
    mutate(cluster_segment = gsub('.*_', '',
                                  cluster_segment)) %>%
    mutate(cluster_segment = case_when(
      grepl('^N[:0-9:]{0,1}$', segment, ignore.case = T) & cluster_segment == 'na' ~ tolower(segment),
      .default = cluster_segment)) %>% 
    filter(tolower(segment) == cluster_segment) %>%
    mutate(cluster_number = case_when(
      is.na(cluster_number) & na.locf0(cluster_number, fromLast = TRUE) == na.locf0(cluster_number, fromLast = FALSE) ~ na.locf0(cluster_number, fromLast = TRUE), 
      .default = cluster_number))
  
  
  return(out)
}
