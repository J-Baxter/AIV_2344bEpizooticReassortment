# Import old meta data
MakeTipNames <- function(data){
  out <- data %>%
    unite(tipnames, 
          virus_subtype,
          clade,
          isolate_id,
          host_order,
          collection_countryname,
          cluster_profile,
          collection_tipdate,
          sep = '|',                        
          remove = FALSE) %>%
    mutate(tipnames = gsub(' ', '_', tipnames))
  
  return(out)
}

RemoveDuplicated <- function(alignment){
  matches <- rownames(alignment) %>% 
    regmatches(., gregexpr("EPI_ISL_\\d+[^.|]*", .)) %>% 
    unlist() 
  
  duplicated_indices <- which(duplicated(matches) | duplicated(matches, fromLast = TRUE))
  duplicated_items <- matches[duplicated_indices]
  duplicated_list <- sapply(duplicated_items, function(x) which(grepl(x, rownames(alignment))), simplify = F)
  
  duplicated_seqnames <- lapply(duplicated_list, function(x) rownames(alignment)[x])
  drop_seqnames <- lapply(duplicated_seqnames, function(x) x[!grepl('^EPI', x)]) %>% unlist()
  
  new_alignment <- alignment[!rownames(alignment) %in% drop_seqnames,]
  return(new_alignment)
}

isDateError <- function(dataframe){
  out <- dataframe %>%
    mutate(collection_dateerror = case_when(
      is.na(collection_tipdate) ~ TRUE,
      as.numeric(gsub('-.*', '', collection_tipdate))>2024 ~ TRUE,
      as.numeric(gsub('-.*', '', collection_tipdate))<1990 ~ TRUE,
      .default = FALSE))
  
  return(out)
}


ReNameAlignment <- function(alignment, data){
  z <- rownames(alignment)
  
  isolates <-  regmatches(z, gregexpr("EPI_ISL_(china_){0,1}\\d+[^.|]*", z)) %>% unlist()
  
  new_seqnames <- sapply(isolates, function(x) data$tipnames[data$isolate_id %in% x]) %>% 
    as.vector() 
  
  rownames(alignment) <-  new_seqnames
  
  return(alignment)
}


#import metadata
metadatafiles <-  list.files('./2023Dec02/metadata',
                             full.names = T, 
                             pattern="csv") %>%
  .[!grepl('n_', .)]


metadata <- lapply(metadatafiles,
                   read_csv,
                   col_types = cols(collection_date = col_date(format = "%Y-%m-%d"),
                                    collection_tipdate = col_character(),
                                    `...28` = col_skip())) 
# Import new alignments
alignmentfiles <- list.files('./2024Jun01/region_alignments_withBLAST',
                             full.names = T, 
                             pattern="fasta")

alignments <- lapply(alignmentfiles,
                     read.dna, 
                     format = 'fasta') %>%
  lapply(., as.matrix) 

no_duplicates_aln <- lapply(alignments, RemoveDuplicated)


# Extract isolate ID

isolate_ids <-  no_duplicates_aln %>%
  lapply(., rownames) %>%
  lapply(., function(x) regmatches(x, gregexpr("EPI_ISL_(china_){0,1}\\d+[^.|]*", x)) %>% unlist()) 


tip_dates <- no_duplicates_aln %>%
  lapply(rownames) %>%
  lapply(function(x) str_extract(x, '(?<=[\\|.])\\d\\d\\d\\d-\\d\\d(-\\d\\d){0,1}'))

# which isolate ID not in meta data
new_isolate_ids <- mapply(function(x,y) y[!y %in% x$isolate_id],
                          metadata,
                          isolate_ids)


# Extract new isolate IDs from new meta data
new_meta <- lapply(new_isolate_ids, function(x) meta %>% filter(isolate_id %in% x))

# format new meta data
temp <- lapply(new_meta, function(x) x %>% as_tibble() %>% unite(., location, starts_with('location'), sep = ' / ', remove = F) %>% 
                               mutate(location = tolower(location)) %>%
                               mutate(location = gsub("[^A-Za-z]", " ", tolower(location))) %>%
                               mutate(location =  gsub("\\s+", " ", str_trim(location))))
  
new_meta_formatted <- list()

for (i in 1:length(temp)){ 
  if(nrow(temp[[i]])>0){
    print(i)
    new_meta_formatted[[i]] <- temp[[i]] %>%
      mutate(across(everything(), .fns = ~ifelse(.x == "", NA, .x))) %>%
      # Format source before taxa allocation
      rename(joint_location = location) %>%
      mutate(across(contains('location'), 
                    .fns = ~ gsub("^NA$|^missing$|^unknown$", NA, .x))) %>%
      mutate(location = coalesce(location_3,location_2)) %>%
      select(-joint_location)%>%
      #rename(source=host) %>%
      # select(-starts_with('host')) %>%
      mutate(date = case_when(is.na(date) & !is.na(date_year) ~ date_year,
                              .default = date)) %>%
      mutate(date_format = case_when(
    grepl('\\d{4}-\\d{2}-\\d{2}', date)  ~ "yyyy-mm-dd", 
    grepl('\\d{2}-\\d{2}-\\d{4}', date) & as.numeric(str_split_i('date', '-', 2)) <= 12 ~ "dd-mm-yyyy", 
    grepl('\\d{4}-\\d{2}', date) ~ "yyyy-mm", 
    grepl('\\d{2}-\\d{4}', date) ~ "mm-yyyy",
    grepl('\\d{4}', date) ~ "yyyy" )) %>%
      
      mutate(date = str_trim(date))%>%
      split(~date_format) %>% 
      map_at("yyyy-mm-dd",  
             ~ mutate(.x,
                      date_parsed = ymd(date),
                      date_ymd = format(date_parsed, '%Y-%m-%d'), 
                      date_dec = decimal_date(date_parsed),
                      date_ym = format(date_parsed, '%Y-%m'),
                      date_y = format(date_parsed, '%Y') )) %>% 
      map_at("dd-mm-yyyy",  
             ~ mutate(.x,
                      date_parsed = dmy(date),
                      date_ymd = format(date_parsed, '%Y-%m-%d'), 
                      date_dec = decimal_date(date_parsed),
                      date_ym = format(date_parsed, '%Y-%m'),
                      date_y = format(date_parsed, '%Y') )) %>% 
      map_at("yyyy-mm",
             ~ mutate(.x,
                      date_parsed = ym(date),
                      date_ymd = NA, 
                      date_dec = NA,
                      date_ym = format(date_parsed,'%Y-%m') ,
                      date_y = format(date_parsed, '%Y') )) %>%
      map_at("mm-yyyy",
             ~ mutate(.x,
                      date_parsed = my(date),
                      date_ymd = NA, 
                      date_dec = NA,
                      date_ym = format(date_parsed,'%Y-%m') ,
                      date_y = format(date_parsed, '%Y') )) %>%
      map_at("yyyy",
             ~ mutate(.x,
                      date_parsed = as.POSIXlt.character(date, format = '%Y'),
                      date_ymd = NA, 
                      date_ym = NA,
                      date_dec = NA,
                      date_y = format(date_parsed, '%Y'))) %>%
      
      list_rbind() %>%
  #rename(month_date = date_year_month) %>%
  #rename(decimal_date = date_frac) %>%
  # mutate(year_date = date_year) %>%
      #select(-matches('^date.+')) %>%
  #rename_with(.fn = ~gsub('_', '.', .x))%>%
  # Location (iso name and subdivision)
      mutate(location = tolower(location)) %>%
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
        grepl('hof am leithaberge|haringsee|fishamend|bruderndorf|tulln|lower austria|nied[(erö)?]sterreich',
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
        grepl("brussels capital region|region de bruxelles[ -]capitale",
              location) ~ "belgium_bruxelles",
        grepl('oost[ -]vlaanderen|west[ -]vlaandere|limburg|vlaams[ -]brabant|antwerp*',
              location) ~ "belgium_vlaanderen",
        grepl("brabant wallon|*hainaut$|*liege$| luxembourg|^namur$|province de namur",
              location) ~ "belgium_wallonie",
        
        # Argentina
        grepl('buenos aires', 
              location) ~ "argentina_buenos aires",
        grepl('catamarca', 
              location) ~ "argentina_catamarca",
        grepl('chaco', 
              location) ~ "argentina_chaco",
        grepl('chubut', 
              location) ~ "argentina_chubut",
        grepl('ciudad de buenos aires',
              location) ~ "argentina_ciudad de buenos aires",
        grepl('c[ó?]rdoba', 
              location) ~ "argentina_córdoba",
        grepl('corrientes',
              location) ~ "argentina_corrientes",
        grepl('entre r[í?]os', 
              location) ~ "argentina_entre ríos",
        grepl('formosa',
              location) ~ "argentina_formosa",
        grepl('jujuy', 
              location) ~ "argentina_jujuy",
        grepl('la pampa', 
              location) ~ "argentina_la pampa",
        grepl('la rioja', 
              location) ~ "argentina_la rioja",
        grepl('mendoza', 
              location) ~ "argentina_mendoza",
        grepl('misiones', 
              location) ~ "argentina_misiones",
        grepl('neuqu[é?]n', 
              location) ~ "argentina_neuquén",
        grepl('r[í?]o negro',
              location) ~ "argentina_río negro",
        grepl('salta', 
              location) ~ "argentina_salta",
        grepl('san juan',
              location) ~ "argentina_san juan",
        grepl('san luis', 
              location) ~ "argentina_san luis",
        grepl('santa cruz',
              location) ~ "argentina_santa cruz",
        grepl('santa fe', 
              location) ~ "argentina_santa fe",
        grepl('santiago del estero', 
              location) ~ "argentina_santiago del estero",
        grepl('tierra del fuego', 
              location) ~ "argentina_tierra del fuego",
        grepl('tucumán', 
              location) ~ "argentina_tucumán",
        
        
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
        grepl('antofagasta', location) ~ "chile_antofagasta",
        grepl('araucan[íi]a', location) ~ "chile_araucanía",
        grepl('arica y parinacota', location) ~ "chile_arica y parinacota",
        grepl('atacama', location) ~ "chile_atacama",
        grepl('aysén del general ibañez del cam', location) ~ "chile_aysén del general ibañez del cam",
        grepl('b[íi]o[- ]{0,1}b[íi]o', location) ~ "chile_bío-bío",
        grepl('coquimbo', location) ~ "chile_coquimbo",
        grepl('libertador general bernardo o\'{0,1}hi', location) ~ "chile_libertador general bernardo o'hi",
        grepl('los lagos', location) ~ "chile_los lagos",
        grepl('los ríos', location) ~ "chile_los ríos",
        grepl('magallanes y antártica chilena', location) ~ "chile_magallanes y antártica chilena",
        grepl('maule', location) ~ "chile_maule",
        grepl('[ñn]uble', location) ~ "chile_ñuble",
        grepl('santiago metropolitan|metropolitana de santiago', location) ~ "chile_santiago metropolitan",
        grepl('tarapac[áa]', location) ~ "chile_tarapacá",
        grepl('valpara[íi]so', location) ~ "chile_valparaíso",
        
        
        
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
        grepl('istarska',
              location) ~ "croatia_istarska",
        grepl('zagreba[cč]ka',
              location) ~ "croatia_zagrebačka", 
        grepl('osjecko-baranjska',
              location) ~ "croatia_osjecko-baranjska", 

        
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
        grepl( "hovedstaden",
              location) ~  "denmark_hovedstaden",
        
        
       
        
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
        grepl('mang(gh){0,1}ystau', 
              location) ~ "kazakhstan_mangghystau",
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
        grepl('tromso|[Tt]roms',
              location) ~"norway_troms", 
        grepl('rogaland',
              location) ~"norway_rogaland",
        grepl("akershus",
               location) ~  "norway_akershus", 
        grepl("ãstfold",
              location) ~  "norway_ãstfold",     
        grepl("aust-agder",
              location) ~  "norway_aust-agder",   
        grepl("buskerud",
              location) ~  "norway_buskerud",      
        grepl("finnmark",
              location) ~  "norway_finnmark",     
        grepl("hedmark",
              location) ~   "norway_hedmark",        
        grepl("hordaland",
              location) ~ "norway_hordaland",  
        grepl("møre og romsdal",
              location) ~ "norway_møre og romsdal", 
        grepl("nord-trøndelag",
              location) ~ "norway_nord-trøndelag",  
        grepl("nordland",
              location) ~ "norway_nordland",        
        grepl("oppland",
              location) ~ "norway_oppland",         
        grepl("oslo",
              location) ~ "norway_oslo",             
        grepl("sogn og fjordane",
              location) ~ "norway_sogn og fjordane", 
        grepl("telemark",
              location) ~ "norway_telemark",          
        grepl("vest-agder",
              location) ~ "norway_vest-agder", 
        grepl("vestfold",
              location) ~  "norway_vestfold",  
        grepl("svalbard",
              location) ~  "svalbard and jan mayen_svalbard", 
        
        
        
        # Peru
        grepl('departamento de lima', 
              location) ~  "peru_lima province",      
        grepl('provincia constitucional del callao',
              location) ~"peru_callao",
        grepl('departamento de tacna',
              location) ~ "peru_tacna",        
        grepl('departamento de ica',
              location) ~  "peru_ica",   
        grepl('lambayeque',
              location) ~  "peru_lambayeque",  
        
        
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
        grepl('jalisco', location) ~ 'méxico_jalisco',
        grepl('guanajuato', location) ~ 'méxico_guanajuato',
        grepl('sonora', location) ~ 'méxico_sonora',
        grepl('veracruz', location) ~ 'méxico_veracruz',
        
        grepl('^choco$', location) ~ 'colombia_chocó',
        grepl('bolivar', location) ~ 'ecuador_bolivar',
        grepl('magdalena', location) ~ 'colombia_magdalena',
        
        
        grepl('esp[íi]rito santo',location) ~ 'brazil_espírito santo',
        grepl('rio de janeiro',location) ~ 'brazil_rio de janeiro',
        grepl('estado do parana|paraná',location) ~ 'brazil_paraná',
       
        
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
        grepl("south moravian region|jihomoravsk[yý] kraj", 
              location) ~"czechia_jihomoravský",
        grepl("kralovehradecky kraj", 
              location) ~"czechia_královéhradecký",
        grepl("liberecky kraj",
              location) ~"czechia_liberecký",
        grepl("mor[oa]vskoslezsk[y?ý] kraj", 
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
        grepl("^aisne$|^nord$|^oise$|pas-de-calais|^somme$|nord[ -]pas[ -]de[ -]calais|picardie", 
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
        
        grepl("bologna|ferrara|forlì-cesena|^modena$|^parma$|piacenza|ravenna|reggio emilia|^rimini$|^fe$|^ra$",
              location) ~ "italy_emilia-romagna",
        
        grepl("gorizia|pordenone|trieste|udine|^ud$", 
              location) ~ "italy_friuli-venezia giulia",
        
        grepl("roma|latina|rieti|viterbo|frosinone",
              location) ~ "italy_lazio",
        
        grepl("genova|imperia|la spezia|savona",
              location) ~ "italy_liguria",
        
        grepl("milano|bergamo|brescia|^bs$|^como$|^cr$|^cremona$|lecco|^lodi$|mantova|^mn$|pavia|sondrio|^varese$|lombardy|^pv$", 
              location) ~ "italy_lombardia",
        
        grepl("ancona|ascoli piceno|fermo|macerata|pesaro e urbino", 
              location) ~ "italy_marche",
        
        grepl("campobasso|isernia", 
              location) ~ "italy_molise",
        
        grepl("torino|alessandria|^asti$|biella|^cuneo$|^novara$|verbano-cusio-ossola|vercelli", 
              location) ~ "italy_piemonte",
        
        grepl("bari|brindisi|foggia|lecce|taranto|barletta-andria-trani|puglia", 
              location) ~ "italy_apulia",
        
        grepl("cagliari|nuoro|oristano|sassari|sud sardegna|sardinia", 
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
        
        grepl("belluno|padova|^pd$|rovigo|treviso|venezia|verona|vicenza|^v[eri]$|^veneto$|venice", 
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
        grepl('^leningrad( oblast){0,1}$|leningrad region', 
              location) ~ "russia_leningrad",
        grepl('^orenburg$', 
              location) ~ "russia_orenburg",
        grepl('^ryazan$', 
              location) ~ "russia_ryazan'",
        grepl('^altai$', 
              location) ~ 'russia_altay',
        
        
        grepl('kherson',
              location) ~ 'ukraine_kherson',
        grepl('p[äa]rnu(maa){0,1}',
              location) ~ 'estonia_pärnu',
        
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
        grepl("england|barrow|northern[_ ]ireland|county of kent|stratford|cheshire|derbyshire|northumberland|lincolnshire|dorset|yorkshire|hertfordshire|gloucestershire|*shire$|norfolk|essex|devon",
              location) ~ "united kingdom",
        grepl("scotland|orkney|western isles|fife|^angus$|dumfries*",
              location) ~ "united kingdom_scotland",
        
        # Japan
        grepl('Kagoshima', 
              isolate_name) ~ 'japan_kagoshima',
        grepl('Chiba', 
              isolate_name) ~ 'japan_chiba',
        grepl('Aichi', 
              isolate_name) ~ 'japan_aichi',
        grepl('Okayama', 
              isolate_name) ~ 'japan_okayama',
        grepl('Niigata',
              isolate_name) ~ 'japan_niigata',
        grepl('Hyogo', 
              isolate_name) ~ 'japan_hyōgo',
        grepl('Ibaraki', 
              isolate_name) ~ 'japan_ibaraki',
        grepl('Kagawa', 
              isolate_name) ~ 'japan_kagawa',
        grepl('Hiroshima', 
              isolate_name) ~ 'japan_hiroshima',
        grepl('Miyazaki', 
              isolate_name) ~ 'japan_miyazaki',
        grepl('Iwate', 
              isolate_name) ~ 'japan_iwate',
        grepl('Hokkaido', 
              isolate_name) ~ 'japan_hokkaido',
       
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
        grepl("catalonia|catalu[nñ]a",
              location) ~ "spain_cataluña",
        

        # Poland
        grepl("lower silesian voivodeship|dolnośląskie", 
              location) ~"poland_dolnośląskie",
        grepl("kuyavian[ -]pomer[ae]nian|kujawsko-pomorskie",
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
        grepl("warmi sko mazurskie|warmian-masurian",
              location) ~"poland_warmińsko-mazurskie",
        grepl("greater poland voivodeship|wielkopolskie", 
              location) ~"poland_wielkopolskie",
        grepl("west pomeranian voivodeship|zachodniopomorskie",
              location) ~"poland_zachodniopomorskie",
        grepl("podkarpackie|subcarpathian voivodeship",
              location) ~"poland_podkarpackie",
        
        
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
        grepl("v[äa]stmanland",
              location) ~"sweden_västmanland",
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
        grepl('^gg$|gyeonggi*|anseong|hwaseong|suwon|gapyeong|yeoncheon|goyang', 
              location) ~ "south korea_gyeonggi-do",
        grepl('^gw$|gangwon-do|cheorwon|goseong', 
              location) ~ "south korea_gangwon-do",
        grepl('^cn$|chungcheongnam*|^asan', 
              location) ~ "south korea_chungcheongnam-do",
        grepl('^cb$|chungcheongbuk*|eumseong|seosan',
              location) ~ "south korea_chungcheongbuk-do",
        grepl('^gn$|gyeongsangnam-do|gimhae',
              location) ~ "south korea_gyeongsangnam-do",
        grepl('^gb$|gyeongsangbuk-do',
              location) ~ "south korea_gyeongsangbuk-do",
        grepl('^jn$|jeollanam-do|jeolla-do|haenam',
              location) ~ "south korea_jeollanam-do",
        grepl('^jb$|jeollabuk-do|iksan|gimje|jeongeup', 
              location) ~ "south korea_jeollabuk-do",
        grepl('^jj$|jeju|seogwipo',
              location) ~ "south korea_jeju",
        grepl('^sw$|seoul', 
              location) ~ "south korea_seoul", 
        grepl('^ic$|incheon',
              location) ~ "south korea_incheon", 
        grepl('^dj$',
              location) ~ "south korea_daejeoni", 
        grepl('^sj$|sejong teugbyeolsi', 
              location) ~ "south korea_sejong", 
        grepl('^gj$|gwangju gwang',
              location) ~ "south korea_gwangju gwang", 
        grepl('^dg$|daegu',
              location) ~ "south korea_daegu", 
        grepl('^ps$|busan', 
              location) ~ "south korea_busan", 
        grepl('^us$|ulsan', 
              location) ~ "south korea_ulsan", 
        
        
        # Netherlands
        grepl('drenthe|westerveld|borger-odoorn|tynaarlo|eext|koekange|middendrexhe|noordenveld|nl assen|meppel',
              location)  ~ 'netherlands_drenthe',
        grepl('flevoland|almere|lelystad|noordoostpolder|noordosterpolder|marker wadden|stroe waddenkust',
              location) ~ 'netherlands_flevoland',
        grepl('zeewolde tulpeiland|zuidlaarmeergebiede oost polder|roggebotsluis|zeewolde|almeerder zand',
              location) ~ 'netherlands_flevoland',
        grepl('friesland|fryslan|tytsjerksteradiel|eastermar|heerenveen|schiermonnikoog|^ameland$|formerum', 
              location) ~ "netherlands_fryslân",
        grepl('kollum|^vlieland.*|ameland|noord[[:space:]]{0,1}buren waddenkust|greonterp|walterswald|nl abbega|terschelling|leeuwarden', 
              location) ~ "netherlands_fryslân",
        grepl('gelderland|zutphen|zevenaar|westvoort|westervoort|andelst|arnhem|bennekom|deliemers',
              location) ~ 'netherlands_gelderland',
        grepl('spijk|doetichem|doetinchem|epe|ermelo|gelmonde|ingen|klarenbeek|lochem|oldenbroek|huissen',
              location) ~ 'netherlands_gelderland',
        grepl('putten|zwolle bergkloosterweg|apeldoorn boogschutterstraat|boven leeuwen|barneveld|wijchen',
              location) ~ 'netherlands_gelderland',
        grepl('groningen|zuidhorn|oldambt|spijk|lauwersmeer|haren boeremapark|scheemda|sint annen|lauwersoog|martunikerk', 
              location)  ~ 'netherlands_groningen',
        grepl('limburg|venlo|gennep|kerkrade|landgraaf|ottersum|maastricht',
              location) ~ 'netherlands_limburg',
        grepl('noord brabant|wernhout|vlijmen|beekendonk|best|boekel|uden|agatha|schaik|rosmalen|^reek&',
              location) ~ 'netherlands_noord-brabant',
        grepl('lierop|oss|north brabant|oirschot|lithse eendenkooi|werkendam|helmond',
              location) ~ 'netherlands_noord-brabant',
        grepl('noord[ -]holland|bloemendaal|heerhugowaard|hilversum|huizen|naarden|reek|texel de schorren|volendam|durgerdam|krommenier', 
              location) ~ 'netherlands_noord-holland',
        grepl('wie[r]{0,1}ingerwerf|den oever|hippolytushoef|kreupel|onderdijk|^obdam$|monnickendam|slootdorp|zuidoost beemster', 
              location) ~ 'netherlands_noord-holland',
        grepl('ijmuiden forteiland|zwaagdijk.*|schellinkhout zuiderdijk|schagen schagerweg|nl wormer$|west graftdijk|oostwoud|hoorn', 
              location) ~ 'netherlands_noord-holland',
        grepl('overijssel|wierden|enschede|hardenberg|heino|losser|overdinkel|raalte|ijsselmuiden|nl zwolle|kamperveen|mastenbroek|markelo',
              location)  ~ 'netherlands_overijssel',
        grepl('utrecht|bilthoven|bosch en duin|bunnik|debilt|derondevenen|soest|haarzuilens|houten|rhenen|eemmeer', 
              location)  ~ 'netherlands_utrecht',
        grepl('ijsselstein|langbroek|nieuwegein|ijsselstein|woerden|leusden|vleuten.*|portengen|eemmeer dode hond|vianen', 
              location)  ~ 'netherlands_utrecht',
        grepl('zeeland|terneuzen|reimerswaal|grenspad|middelburg|weversinlaag|weversinlaag|neeltje jans', 
              location) ~ 'netherlands_zeeland',
        grepl('walcheren t vroon westkapelle|sloehaven', 
              location) ~ 'netherlands_zeeland',
        grepl('zuid holland|zoeterwoude|westland|den haag|rotterdam|rijswijk|leiden|oud albas|oud alblas',
              location) ~ 'netherlands_zuid-holland',
        grepl('nederlek|hardinxveldgiessendam|south holland|eendenkooi stolwijk|haringvliet|bliek$|maasdam', 
              location) ~ 'netherlands_zuid-holland',
        grepl('kwade hoek stellendam|stellendam scheelhoekeiland|leidschendam|nl gouda|groot ammers|reeuwijk|stolwijk', 
              location) ~ 'netherlands_zuid-holland',
        grepl('^nl$|nl sovon|eendenkooi slijkerman', location) ~ 'netherlands',
        
        
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
        grepl('republika srpska', 
              location)~ "bosnia and herzegovina_repuplika srpska",
        
        grepl('sva |^sva$|north america|zeebrugge belgi',
              location) ~ NA,
        
        grepl('central luzon', 
              location) ~"philippines_nueva ecija",
        
        grepl('kokkola',
              location) ~  'finland_western finland',
       
        # Indonesia
        grepl("hulu sungai utara|banjarbaru|south kalimantan", 
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
        
        
        .default = location)) %>%
      left_join(geodata, by = join_by(location == match)) %>% 
      as_tibble() %>%
      
      
      # Format column names
      #dplyr::select(-c(virus_species, id_unsure)) %>%
      dplyr::rename(
        virus_subtype = subtype,
        collection_date = date_ymd,
        collection_datedecimal = date_dec,
        collection_datemonth = date_ym,
        collection_dateyear = date_y,
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
        #host_class = class,
        #host_order = order,
        #host_family = family,
        #host_sciname = sci_name,
        #host_commonname = primary_com_name,
        #host_isbird = is_bird,
        host_isdomestic = domestic_wild
      ) %>%
      
      # Infer simplified host categories for BEAST
      mutate(across(starts_with('host'), .fns = ~tolower(.x))) %>%
      mutate(host_simplifiedhost = case_when(
        host_order %in% c('anseriformes', 'galliformes', 'charadriiformes') ~ paste(host_order, 
                                                                                    host_isdomestic,
                                                                                    sep = '-'),
        host_class == 'mammalia' ~ 'mammal',
        host_class == 'mammal' ~ 'mammal',
        #host_commonname == 'unknown' ~ 'unknown',
        #host_commonname == 'environment' ~ 'environment',
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
        #host_family,
        # host_sciname,
        #host_commonname,
        #host_isbird,
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
      #select(-c(location, varname_1, iso_1)) %>%
      
      # Format NAs
      mutate(across(everything(), .fns = ~ gsub('^NA$', NA, .x))) 
    
    
   
    
    
  }else{
    new_meta_formatted[[i]] <-tibble(virus_subtype = character(),
                                     isolate_id = character(),
                                     isolate_name = character(),
                                     collection_date = character(),
                                     collection_datedecimal = character(),
                                     collection_datemonth = character(),
                                     collection_dateyear = character(),
                                     host_class = character(),
                                     host_order = character(),
                                     host_isdomestic = character(),
                                     host_simplifiedhost = character(),
                                     collection_regionname = character(),
                                     collection_countryname = character(),
                                     collection_countrycode = character(),
                                     collection_countrylat = character(),
                                     collection_countrylong = character(),
                                     collection_subdiv1name = character(),
                                     collection_subdiv1code = character(),
                                     collection_subdiv1lat = character(),
                                     collection_subdiv1long = character(),
                                     clade = character(),
                                     location_1 = character(),
                                     location_2 = character(),
                                     location_3 = character(),
                                     location_4 = character(),
                                     location_5 = character(),
                                     submission_date = character(),
                                     date = character(),
                                     date_frac = character(),
                                     date_year = character(),
                                     cluster_pb2 = character(),
                                     cluster_pb1 = character(),
                                     cluster_pa = character(),
                                     cluster_ha = character(),
                                     cluster_np = character(),
                                     cluster_na = character(),
                                     cluster_mp = character(),
                                     cluster_ns = character(),
                                     cluster_profile = character(),
                                     cluster_label = character(),
                                     location = character(),
                                     date_format = character(),
                                     date_parsed = character(),
                                     ...1 = character(),
                                     varname_1 = character(),
                                     iso_1 = character()
                                     )
    
      
    
  }
  
}


#new_locations <- test_all %>% select(contains('location')) %>% distinct() 
stopifnot(all( lapply(new_meta_formatted, nrow) %>%
                 unlist() - lapply(new_isolate_ids, length) %>% unlist() == 0))


new_meta_formatted <- lapply(new_meta_formatted, function(x) if(!is.null(x)){x %>% 
                               select(-c(starts_with('location'),
                                         cluster_pb1,
                                         cluster_pb2,
                                         cluster_pa,
                                         cluster_ha,
                                         cluster_np,
                                         cluster_na,
                                         cluster_mp,
                                         cluster_ns,
                                         varname_1,
                                         iso_1,
                                         "...1"))})

old_isolate_ids <- mapply(function(x,y) x %>% filter(isolate_id %in% y),
                          metadata,
                          isolate_ids,
                          SIMPLIFY = F) %>% 
  
  lapply(., function(x) x[!!rowSums(!is.na(x)),])

        

# combine all metadata

combined_metadata <- mapply(function(x,y) list(x,y) %>% 
                              lapply(., function(z) z %>% 
                                       mutate(across(everything(), .fns = ~ as.character(.x)))) %>% 
                              bind_rows(.) ,
                            old_isolate_ids,
                            new_meta_formatted,
                            SIMPLIFY = F) %>%
  lapply(., function(x) x %>%
           select(-c(genome,
                     cluster_profileclass,
                     cluster_profilecolour, 
                     collection_tipdate,
                     segment,
                     collection_dateerror,
                     cluster_segment, 
                     cluster_number,
                     submission_date,
                     date,date_frac,
                     date_year, 
                     date_format,
                     date_parsed, 
                     virus_species,
                     week_date)) %>%
           mutate(cluster_profile = coalesce(profile, cluster_profile)) %>%  
           mutate(cluster_label = coalesce(profile_lab, cluster_label)) %>%
           select(-starts_with('profile'))) %>%
  lapply(., function(x) x %>% 
           mutate(collection_tipdate = coalesce(collection_date,collection_datemonth, collection_dateyear))) %>%
  lapply(.,function(x) x %>%
           mutate(clade = gsub('[[:punct:]]', '', clade))) %>%
  lapply(., MakeTipNames)
  

# sanity check nrow(meta) == nrow(aln) & all(isolate_id %in% "isolate id's from name")
all(lapply(combined_metadata, nrow) %>% unlist() == lapply(no_duplicates_aln , nrow) %>% unlist())


# assign new seqnames to alignment
renamed_aln <- mapply(ReNameAlignment,
                      no_duplicates_aln,
                      combined_metadata,
                      SIMPLIFY = F)


# Segment names
segnames <- str_split(alignmentfiles,  '_') %>% 
  lapply(., tail, n = 3) %>% 
  lapply(., `[`, 2) %>%
  lapply(., function(x) gsub("\\..*$", "",x)) %>%
  lapply(., paste0 , collapse = '_') %>%
  unlist() %>%
  tolower() %>%
  gsub('na', 'n', .)


combined_metadata_imp <- mapply(function(x,y) x %>% mutate(segment = y), combined_metadata, as.list(segnames), SIMPLIFY = F) %>%
  mapply(ImputeCladeandCluster,  ., renamed_aln, SIMPLIFY = F)   %>%
  lapply(., function(x) x %>% mutate(cluster_number = paste0('profile', str_pad(cluster_number, 3, pad = "0"))))
  


####################################################################################################
# Check year and remove if problematic
combined_metadata_imp_dateschecked <- combined_metadata_imp %>%
  lapply(., function(x) x %>% 
           isDateError() %>%
           filter(!collection_dateerror))

prob_seqs <- lapply(combined_metadata_imp, function(x) x %>% 
                      isDateError() %>%
                      filter(collection_dateerror) %>% 
                      pull(tipnames))

renamed_aln_dateschecked <- renamed_aln %>%
  mapply(function(alignment,seqnames_to_drop) alignment[!rownames(alignment) %in% seqnames_to_drop,], 
         .,
         prob_seqs, 
         SIMPLIFY = FALSE)
####################################################################################################
# Alignments and meta data to file
alignmentfiles_newnames <- alignmentfiles  %>%
  gsub('./2024Jun01/region_alignments_withBLAST/', './2024Jul12/region_alignments/',. ) 


mapply(write.dna, 
       renamed_aln_dateschecked , 
       alignmentfiles_newnames,
       format = 'fasta')


alignmentfiles_newnames <- metadatafiles %>%
  gsub('./2024Jun01/metadata/', './2024Jul12/region_metadata/',. ) 

mapply(write_csv, 
       quote= 'needed',
       combined_metadata_imp_dateschecked, 
       metadatafiles)  



