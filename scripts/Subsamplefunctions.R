# Data Import and Subsampling Functions
FormatTipData <- function(x) {
  isolate.id <- x[grep("EPI_ISL_*", x)]
  
  subtype <- x[grep("H[[:digit:]]{1,}N[[:digit:]]*", x)] %>%
    str_extract(., "H[[:digit:]]{1,}N[[:digit:]]*")
  
  date <- x[grep("^[[:digit:]][[:digit:]][[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]-[[:digit:]][[:digit:]]", x)]
  
  isolate.name <- x[grep("/[^/]*/[^/]*/", x)]
  
  if (identical(isolate.id, character(0))) {
    isolate.id <- NA
  }
  
  if (identical(subtype, character(0))) {
    subtype <- NA
  }
  
  if (identical(date, character(0))) {
    date <- NA
  }
  
  if (identical(isolate.name, character(0))) {
    isolate.name <- NA
  }
  
  out <- cbind.data.frame(
    "isolate.id" = isolate.id,
    "subtype" = subtype,
    "date" = date,
    "isolate.name" = isolate.name
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
    mutate(decimal.date = format(round(decimal_date(date), 2), 
                                 nsmall = 2) ) %>%
    mutate(decimal.date = suppressWarnings(as.double(decimal.date))) %>%
    mutate(week.date = format(date, "%Y-%V")) %>%
    mutate(month.date = format(date, "%Y-%m")) %>%
    
    # Format isolate name
    mutate(isolate.name = gsub("_(\\d{4})|(\\d{4})_", "\\1", perl = TRUE, isolate.name)) %>%
    mutate(isolate.name = case_when(!grepl("/(\\d{4})$", isolate.name) & !grepl("/$", isolate.name) ~ 
                                      paste0(isolate.name, "/NA"),
                                    !grepl("/(\\d{4})$", isolate.name) & grepl("/$", isolate.name) ~ paste0(isolate.name, "NA"),
                                    .default = isolate.name
    )) %>%
    
    # Extract data from isolate name
    separate_wider_delim(isolate.name,
                         delim = "/",
                         names = c("virus_species", "source", "location", "id_unsure", "year"),
                         cols_remove = F,
                         too_few = "align_end",
                         too_many = "debug"
    ) %>%
    
    # Drop extraneous cols and format
    dplyr::select(-c(year, contains("isolate.name_"))) %>%
    mutate(virus_species = case_when(str_length(source) == 1 ~ source,
                                     .default = virus_species
    )) %>%
    mutate(source = case_when(str_length(source) == 1 ~ NA,
                              .default = source
    )) 
  
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
      rename(month.date = date_year_month) %>%
      rename(decimal.date = date_frac)
    
  } else{
    data <- data %>%
      mutate(date = parsedate::parse_date(date) %>%
               as.Date()) %>%
      mutate(decimal.date = format(round(decimal_date(date), 2), 
                                   nsmall = 2) ) %>%
      mutate(decimal.date = suppressWarnings(as.double(decimal.date))) %>%
      mutate(month.date = format(date, "%Y-%m"))
  }
  
  data <- data %>%
    select(-matches('^date.+')) %>%
    rename_with(.fn = ~gsub('_', '.', .x))%>%
    # Location (iso name and subdivision)
    mutate(location = tolower(location)) %>%
    mutate(location = gsub("[^A-Za-z]", " ", tolower(location))) %>%
    mutate(location =  gsub("\\s+", " ", str_trim(location))) %>%
    mutate(location = case_when(
      
      # Austria
      grepl('burgenland|stegersbach|mattersburg|bad sauerbrunn',
            location) ~ 'burgenland',
      grepl('karnten', 
            location)  ~ 'karnten',
      grepl('niederosterreich|voesendorf|strasshof|stockerau|sollenau|breitensee|biberbach',
            location) ~ "niederösterreich", 
      grepl('sankt poelten|pottschach|neunkirchen|markgrafneusiedl|koenigstetten|klosterneuburg',
            location) ~ "niederösterreich",
      grepl('hof am leithaberge|haringsee|fishamend|bruderndorf|tulln|lower austria',
            location) ~ "niederösterreich",
      grepl('oberosterreich$|sankt georgen|^linz$|oberneukirchen',
            location) ~ "oberösterreich",
      grepl('salzburg', 
            location)  ~ 'salzburg',
      grepl('steiermark|tagnitz|^graz$|eggersdorf|leibnitz', 
            location) ~ "steiermark",
      grepl('tirol$', 
            location)  ~ 'tirol',
      grepl('vorarlberg|voralrberg', 
            location)  ~ "vorarlberg",
      grepl('^wien$|^vienna$|bundesland wien', 
            location) ~ 'wien',
      grepl('carinthia', 
            location) ~ "kärnten",
      
      # Belgium
      grepl('antwerpen',
            location) ~ "antwerpen",
      grepl("brabant wallon",
            location) ~ "brabant wallon",
      grepl("brussels capital region|region de bruxelles capitale",
            location) ~ "brussels-capital region",
      grepl('*hainaut$',
            location) ~ "hainaut",
      grepl('*liege$',
            location) ~ "liege",
      grepl('limburg',
            location) ~ "limburg",
      grepl(' luxembourg', 
            location) ~ "luxembourg",
      grepl("^namur$|province de namur",
            location) ~ "namur",
      grepl('oost vlaanderen',
            location) ~ "oost-vlaanderen",
      grepl('vlaams brabant', 
            location) ~ "vlaams-brabant",
      grepl('west vlaanderen', 
            location) ~ "west-vlaanderen",
      
      # Bulgaria 
      grepl('blagoevgrad', 
            location) ~ "blagoevgrad",
      grepl('burgas$', 
            location) ~ "burgas",
      grepl('dobrich$', 
            location) ~ "dobrich",
      grepl('gabrovo$', 
            location) ~ "gabrovo",
      grepl('haskovo$', 
            location) ~ "haskovo",
      grepl('kardzhali$', 
            location) ~ "kardzhali",
      grepl('kyustendil$', 
            location) ~ "kyustendil",
      grepl('lovech$', 
            location) ~ "lovech",
      grepl('blagoevgrad', 
            location) ~ "blagoevgrad",
      grepl('"pazardzhik"', 
            location) ~ "pazardzhik",
      grepl('pernik', 
            location) ~ "pernik",
      grepl('pleven|slavyanovo', 
            location) ~ "pleven",
      grepl('plovdiv$', 
            location) ~ "plovdiv",
      grepl('^razgrad$', 
            location) ~ "razgrad",
      grepl('^ruse$', 
            location) ~ "ruse",
      grepl('shumen', 
            location) ~ "shumen",
      grepl('silistra', 
            location) ~ "silistra",
      grepl('sliven', 
            location) ~ "sliven",
      grepl('smolyan', 
            location) ~ "smolyan",
      grepl('^sofia$', 
            location) ~ "sofia",
      grepl('sofia-grad', 
            location) ~ "sofia-grad",
      grepl('stara zagora', 
            location) ~ "stara zagora",
      grepl('targovishte', 
            location) ~ "targovishte",
      grepl('varna', 
            location) ~ "varna",
      grepl('veliko tarnovo', 
            location) ~ "veliko tarnovo",
      grepl('^vidin$', 
            location) ~ "vidin",
      grepl('^vratsa$', 
            location) ~ "vratsa",
      grepl('^yambol$', 
            location) ~ "yambol",
      
      # Burkina Faso
      grepl("koubri", 
            location) ~ "kadiogo",
      grepl("gomboussougou", 
            location) ~ "zoundwéogo",
      
      # Cambodia
      grepl("krong phnum penh", 
            location) ~"phnom penh",
      grepl("khett kandal",
            location) ~"kandal",
      grepl("khett prey veng",
            location) ~"prey veaeng",
      grepl("takeo", 
            location) ~"taakaev",
      
      # Chile
      grepl("region de antofagasta", 
            location) ~ "antofagasta",
      
      # Colombia
      grepl("departamento de bolivar", 
            location) ~"bolívar",
      grepl("departamento del choco", 
            location) ~"chocó",
      grepl("departamento de cordoba",
            location) ~"córdoba",
      grepl("departamento del magdalena", 
            location) ~"magdalena",
      
      # Costa Rica
      grepl("provincia de limon",
            location) ~ "limón",
      
      # Croatia
      grepl('koprivnicko krizevacka zupanija', 
            location) ~ "koprivničko-križevačka županija",
      grepl('vukovarsko srijemska zupanija',
            location) ~ "vukovarsko-srijemska županija",
      grepl('sisacko moslavacka zupanija',
            location) ~ "sisačko-moslavačka županija",   
      
      # Denmark
      grepl('region sjaland', 
            location) ~    "sjælland",
      grepl('region syddanmark', 
            location) ~ "syddanmark",
      grepl('region nordjylland', 
            location) ~ "nordjylland",
      grepl('region midtjylland',
            location) ~"midtjylland",
      
      # Ecuador
      grepl('provincia de manabi', 
            location) ~ "manabí",
      
      # Honduras
      grepl('departamento de atlantida',
            location) ~ "atlántida",
      
      # Hungary
      grepl('bacs kiskun',
            location) ~   "bács-kiskun",
      grepl('csongrad megye', 
            location) ~ "csongrád",
      grepl('bekes megye', 
            location) ~ "békés",
      grepl('hajdu bihar', 
            location) ~ "hajdú-bihar",
      grepl('somogy megye', 
            location) ~  "somogy",
      
      # Iceland
      grepl('westfjords',
            location) ~ "vestfirðir",
      
      # Indonesia
      grepl('east java', 
            location) ~ "jawa timur",
      
      # Iraq
      grepl('muhafazat ninawa',
            location) ~ "ninawa",
      
      # Israel
      grepl('northern district', 
            location) ~ "hazafon",
      
      # Kazakhstan
      grepl('akmola', 
            location) ~ "aqmola oblysy", 
      grepl('kostanay', 
            location) ~ "qostanay oblysy",        
      grepl('north kazakhstan',
            location) ~ "soltüstik quzaqstan oblysy",
      grepl('almaty province', 
            location) ~ "almaty oblysy",
      
      # Latvia
      grepl('jurmala',
            location) ~"jūrmala",
      # Lithuania
      grepl('vilnius', 
            location) ~"vilniaus apskritis",
      # Luxembourg
      grepl('district de grevenmacher', 
            location) ~ "grevenmacher",
      grepl('district de diekirch', 
            location) ~ "diekirch",
      # Mali
      grepl('kati', 
            location) ~ "koulikoro",
      
      # Montenegro
      grepl('skadar lake',
            location) ~ "bar",
      
      # Niger
      grepl('bouza tahoua', 
            location) ~ "tahoua",
      grepl('torodi tillab ri',
            location) ~ "tillabéri",  
      grepl(' zinder$', 
            location) ~ "zinder",
      grepl('maradi$',
            location) ~ "maradi",
      grepl('niamey$', 
            location) ~ "niamey",
      
      # Norway
      grepl('hitra', 
            location) ~ "sør-trøndelag",
      grepl('tromso',
            location) ~"troms",
      
      # Peru
      grepl('departamento de lima', 
            location) ~  "municipalidad metropolitana de lima",      
      grepl('provincia constitucional del callao',
            location) ~"el callao",
      grepl('departamento de tacna',
            location) ~ "tacna",        
      grepl('departamento de ica',
            location) ~  "ica",             
      
      # Portugal
      grepl('distrito de setubal', 
            location) ~ "setúbal",
      
      # Romania
      grepl('judetul timis', 
            location) ~ "timiș",     
      grepl('judetul ilfov', 
            location) ~ "ilfov",      
      grepl('judetul constanta', 
            location) ~ "constanța", 
      grepl('mures ungheni', 
            location) ~ "mureș",     
      grepl('judetul giurgiu', 
            location) ~ "giurgiu",   
      grepl('tulcea',
            location) ~ "tulcea",    
      
      # Senegal
      grepl('region de thies',
            location) ~ 'thiès',
      grepl('djoudj bird sanctuary',
            location) ~ 'saint-louis',
      
      # Slovakia
      grepl('dobrohost', 
            location) ~ "trnavský kraj",
      grepl('kosice',
            location) ~  "košický kraj" ,
      grepl('kalinkovo|petrzalka',
            location) ~"bratislavský kraj",
      
      # South Africa
      grepl('province of the western cape',
            location) ~ "western cape",
      grepl('province of kwazulu natal', 
            location) ~  "kwazulu-natal",
      grepl('province of north west',
            location) ~ "north-west (south africa)",
      
      # Switzerland
      grepl('kanton zurich',
            location) ~ "zürich",
      grepl('cantone ticino', 
            location) ~ "ticino",    
      grepl('kanton basel stadt', 
            location) ~ "basel-stadt",
      grepl('kanton luzern', 
            location) ~ "luzern",
      
      # Venezuela
      grepl("estado anzoategui", 
            location) ~ "anzoátegui",
      
      
      # Vietnam
      grepl("vietnam", 
            location) ~ "viet nam",
      grepl("quang ninh" ,
            location) ~  "quảng ninh" ,
      grepl('tinh lang son',
            location) ~ "lạng sơn",
      grepl('tinh thanh hoa',
            location) ~ "thanh hóa",
      grepl('nghe an|nhge an', 
            location) ~  "nghệ an", 
      grepl('tinh quang tri',
            location) ~ "quảng trị",
      grepl('nam ha tinh', 
            location) ~ "hà tỉnh",
      grepl('thu do ha noi',
            location) ~ "hà nội, thủ đô",
      
      
      # PRC
      grepl("northern china|peking", 
            location) ~ "beijing", # Interpretation from paper - check with Lu
      grepl("chongqin*", 
            location) ~ "chongqing",
      grepl("shanghai",
            location) ~ "shanghai",
      grepl('tianjing',
            location) ~ 'tianjin',
      grepl('anhui',
            location) ~ 'anhui',
      grepl('^fujain$|^fujian$', 
            location) ~ 'fujian',
      grepl('gansu',
            location) ~ 'gansu',
      grepl('qingyuan|foshan|guangdong',
            location) ~ 'guangdong',
      grepl('guizhou', 
            location) ~ 'guizhou',
      grepl("heinan|hainan", 
            location) ~ "hainan", # confirm 
      grepl("wuhan|hebei", 
            location) ~ "hebei",
      grepl("heilongjiang", 
            location) ~ "heilongjiang",
      grepl("henan|sanmenxia|dongting", 
            location) ~ "henan",
      grepl("hubei",
            location) ~ "hubei",
      grepl("hunan|changsha", 
            location) ~ "hunan",
      grepl("jiangsu|suzhou|xuzhou",
            location) ~ "jiangsu",
      grepl("jiangxi", 
            location) ~ "jiangxi",
      grepl("jilin", 
            location) ~ "jilin",
      grepl("liaoning",
            location) ~ "liaoning",
      grepl("qinghai", 
            location) ~ "qinghai",
      grepl("shaanxi", 
            location) ~ "shaanxi",
      grepl("shandong|rongcheng", 
            location) ~ "shandong",
      grepl("shanxi|southwestern china",
            location) ~ "shanxi",
      grepl("sichuan",
            location) ~ "sichuan",
      grepl("taiwan", 
            location) ~ "taiwan",
      grepl("yunnan.*", 
            location) ~ "yunnan",
      grepl("zhejiang|hangzhou|nanji$|eastern china",
            location) ~ "zhejiang",
      grepl("guangxi$|guilin|hechi|guangxi zhuang",
            location) ~ "guangxi",
      grepl("inner mongolia|tumuji|nei mongol", 
            location) ~ "nei mongol",
      grepl("ningxia$|ningxia hui", 
            location) ~ "ningxia",
      grepl("xinjiang$", 
            location) ~ "xinjiang",
      grepl("xizang$|tibet", 
            location) ~ "xizang",
      grepl("macau",
            location) ~ "aomen (macau)",
      grepl("hong kong|hongkong",
            location) ~ "xianggang (hong-kong)",
      
      
      # Canada
      grepl("^ab$",
            location) ~ "alberta",
      grepl("^bc$", 
            location) ~ "british columbia",
      grepl('^mb$', 
            location) ~'manitoba',
      
      # Czechia
      grepl("jihocesky kraj", 
            location) ~"jihočeský kraj",
      grepl("south moravian region", 
            location) ~"jihomoravský kraj",
      grepl("kralovehradecky kraj", 
            location) ~"královéhradecký kraj",
      grepl("liberecky kraj",
            location) ~"liberecký kraj",
      grepl("morovskoslezsky kraj", 
            location) ~"moravskoslezský kraj",
      grepl("olomoucky kraj",
            location) ~"olomoucký kraj",
      grepl("pardubicky kraj", 
            location) ~"pardubický kraj",
      grepl("plzensky kraj",
            location) ~"plzeňský kraj",
      grepl("hlavni mesto praha", 
            location) ~"praha, hlavní město",
      grepl("stredocesky kraj",
            location) ~"středočeský kraj",
      grepl("ustecky kraj",
            location) ~"ústecký kraj",
      grepl("vysocina kraj",
            location) ~"vysočina",
      grepl("zlinsky kraj",
            location) ~"zlínský kraj",
      grepl("nov hrady by ov", 
            location) ~"české budějovice",
      
      
      # France
      grepl("^france&",
            location) ~"france",
      grepl("^alsace&",
            location) ~"alsace",
      grepl("aquitaine",
            location) ~"aquitaine",
      grepl("auvergne",
            location) ~"auvergne",
      grepl("basse normandie",
            location) ~"basse-normandie",
      grepl("bourgogne", 
            location) ~"bourgogne",
      grepl("brittany",
            location) ~"bretagne",
      grepl("^centre$",
            location) ~"centre",
      grepl("champagne ardenne",
            location) ~"champagne-ardenne",
      grepl("corsica", 
            location) ~"corse",
      grepl("franche comté",
            location) ~"franche-comté",
      grepl("haute normandie", 
            location) ~"haute-normandie",
      grepl("region ile de france",
            location) ~"île-de-france",
      grepl("languedoc roussillon", 
            location) ~"languedoc-roussillon",
      grepl("limousin", 
            location) ~"limousin",
      grepl("region lorraine", 
            location) ~"lorraine",
      grepl("midi pyrenees", 
            location) ~"midi-pyrénées",
      grepl("region nord pas de calais", 
            location) ~"nord - pas-de-calais",
      grepl("region pays de la loire", 
            location) ~"pays de la loire",
      grepl("poitou charentes",
            location) ~"poitou-charentes",
      grepl("rhone alpes", 
            location) ~"rhône-alpes",
      grepl("guadeloupe", 
            location) ~"guadeloupe",
      grepl("guyane", 
            location) ~"guyane",
      grepl("martinique", 
            location) ~"martinique",
      grepl("réunion",
            location) ~"réunion",
      grepl("ain", 
            location) ~"ain",
      grepl("aisne", 
            location) ~"aisne",
      grepl("allier",
            location) ~"allier",
      
      
      # Italy
      grepl("emilia romagna", 
            location) ~"emilia-romagna",
      grepl("lombardy",
            location) ~"lombardia",
      grepl("^bs$", 
            location) ~"brescia",
      grepl("^cr$", 
            location) ~"cremona",
      grepl("^mn$", 
            location) ~"mantova",
      grepl("^pd$",
            location) ~"padova",
      grepl("^ud$",
            location) ~"udine",
      grepl("^ve$", 
            location) ~"venezia",
      grepl("^vr$", 
            location) ~"verona",
      grepl("^vi$", 
            location) ~"vicenza",
      
      # Kosovo (Currently represented as serbia)
      grepl('gjilan',
            location) ~ "kosovsko-pomoravski okrug",  
      grepl('mitrovice|mitrovica', 
            location) ~ "kosovsko-mitrovački okrug",
      grepl('shtime|rahovec|malisheve',
            location) ~ "pećki okrug",
      grepl('kacanik|ferizaj|shtime|suhareke',
            location) ~ "kosovski okrug",
      grepl('kosovo', 
            location) ~   "kosovo-metohija",
      
      # Russian Federation
      grepl('russian federation', 
            location) ~ 'russia',
      grepl('khabarovsk', 
            location) ~ "khabarovskiy kray",
      grepl('chelyabinsk', 
            location) ~ "chelyabinskaya oblast'",
      grepl("astrakhan", 
            location) ~ "astrakhanskaya oblast'",
      grepl("novosibirsk|chany lake|^chany$", 
            location) ~ "novosibirskaya oblast'",
      grepl("omsk.*", 
            location) ~ "omskaya oblast'",
      grepl("rostov oblast|rostov on don", 
            location) ~ "rostovskaya oblast'",
      grepl("russia primorje", 
            location) ~ "primorskiy kray",
      grepl("sakhalin", 
            location) ~ "sakhalinskaya oblast'",
      grepl("yakutia", 
            location) ~ "sakha, respublika [yakutiya]",
      grepl("kostroma", 
            location) ~ "kostromskaya oblast'",
      grepl("amur region", 
            location) ~ "amurskaya oblast'",
      grepl("buryatia", 
            location) ~ "buryatiya, respublika",
      grepl('north ossetia-alania', 
            location) ~ 'severnaya osetiya-alaniya, respublika',
      grepl('dagestan', 
            location) ~ 'dagestan, respublika',
      grepl('stavropol', 
            location) ~ "stavropol'skiy kray",
      grepl('krasnodar', 
            location) ~ "krasnodarskiy kray", 
      grepl('tyumen', 
            location) ~ "tyumenskaya oblast'",
      grepl("magadan",
            location) ~ "magadanskaya oblast'",
      grepl('saratov', 
            location) ~ "saratovskaya oblast'",
      grepl('kurgan*', 
            location) ~ "kurganskaya oblast'",
      grepl('central russia',
            location) ~ 'russia',
      grepl('republic of tatarstan', 
            location) ~"tatarstan, respublika",
      grepl('republic of kalmykia', 
            location) ~ "kalmykiya, respublika",
      
      # USA
      grepl("north dakota",
            location) ~ "north dakota",
      grepl("massachussetts",
            location) ~"connecticut",
      grepl("conneticut",
            location) ~"massachusetts",
      grepl("deleware", 
            location) ~"delaware",
      
      
      # UK
      grepl("england|wales",
            location) ~ "england and wales",
      grepl("county of kent", 
            location) ~ "kent",
      grepl("county borough of wrexham", 
            location) ~ "wrexham;wrecsam",
      grepl("stratford", 
            location) ~ "warwickshire",
      grepl("anglesey",
            location) ~ "isle of anglesey;sir ynys môn",
      grepl("orkney", 
            location) ~ "orkney islands",
      grepl("western isles",
            location) ~ "highland",
      grepl("cheshire", 
            location) ~ "cheshire east",
      grepl("scotland",
            location) ~ "scotland",
      
      # Japan
      grepl("tsukuba", 
            location) ~ "ibaraki",
      
      # Spain
      grepl("castilla[ ]{0,1}la[ ]{0,1}mancha|castille[ ]{0,1}la[ ]{0,1}mancha",
            location) ~ "castilla-la mancha",
      
      # Poland
      grepl("lower silesian voivodeship", 
            location) ~"dolnośląskie",
      grepl("kuyavian pomeranian voivodeship",
            location) ~"kujawsko-pomorskie",
      grepl("lublin voivodeship", 
            location) ~"lubelskie",
      grepl("lubusz voivodeship",
            location) ~"lubuskie",
      grepl("lodzkie|odz voivodeship",
            location) ~"łódzkie",
      grepl("lesser poland voivodeship", 
            location) ~"małopolskie",
      grepl("masovian voivodeship",
            location) ~"mazowieckie",
      grepl("opole voivodeship", 
            location) ~"opolskie",
      grepl("pomeranian voivodeship",
            location) ~"pomorskie",
      grepl("silesian voivodeship", 
            location) ~"śląskie",
      grepl("swietokrzyskie voivodeship", 
            location) ~"świętokrzyskie",
      grepl("warmi sko mazurskie",
            location) ~"warmińsko-mazurskie",
      grepl("greater poland voivodeship", 
            location) ~"wielkopolskie",
      grepl("west pomeranian voivodeship",
            location) ~"zachodniopomorskie",
      
      
      # Sweden
      grepl("blekinge lan",
            location) ~"blekinge län",
      grepl("dalarnas lan",
            location) ~"dalarnas län",
      grepl("gotlands lan",
            location) ~"gotlands län",
      grepl("hallands lan", 
            location) ~"hallands län",
      grepl("jonkopings lan"
            , location) ~"jönköpings län",
      grepl("kalmar lan",
            location) ~"kalmar län",
      grepl("skane lan",
            location) ~"skåne län",
      grepl("stockholms lan",
            location) ~"stockholms län",
      grepl("sodermanlands lan", 
            location) ~"södermanlands län",
      grepl("uppsala lan", 
            location) ~"uppsala län",
      grepl("vastra gotalands lan",
            location) ~"västra götalands län",
      grepl("ostergotlands lan", 
            location) ~"östergötlands län",
      
      
      # Germany 
      grepl("^germany$", 
            location) ~"germany",
      grepl("germany[- ]bw|baden wuerttemberg|freiburg im breisgau|rosengarten|stutensee",
            location) ~"baden-württemberg",
      grepl("germany[- ]by|bavaria|munich|nuremberg|wuerzburg", 
            location) ~"bayern",
      grepl("germany[- ]hb", 
            location) ~"bremen",
      grepl("germany[- ]hh", 
            location) ~"hamburg",
      grepl("germany[- ]he|hessen|freigericht|giessen|karben|ranstadt", 
            location) ~"hessen",
      grepl("germany[- ]ni|lower saxony|brietlingen|dinklage|hannover|lueneburg|nordhorn|osnabruck",
            location) ~"niedersachsen", 
      grepl("osterholz[- ]scharnbeck|theene|wingst",
            location) ~"niedersachsen",
      grepl("germany[- ]nw|north rhine[- ]westphalia|aachen|brueggen", 
            location) ~"nordrhein-westfalen",
      grepl("germany[- ]rp|rhineland[- ]palatinate", 
            location) ~"rheinland-pfalz",
      grepl("germany[- ]sl", 
            location) ~"saarland",
      grepl("germany[- ]sh|schleswig[- ]holstein|aumuehle|brickeln|luebeck|pinneberg|prohn", 
            location) ~"schleswig-holstein",
      grepl("germany[- ]be",
            location) ~"berlin",
      grepl("germany[- ]bb", 
            location) ~"brandenburg",
      grepl("germany[- ]mv|mecklenburg[- ]vorpommern|grevesmuehlen|parchim|ruegen|warin", 
            location) ~"mecklenburg-vorpommern",
      grepl("germany[- ]sn|sachsen|doberschuetz|dresden|gross dueben|leipzig|saxony",
            location) ~"sachsen",
      grepl("germany[- ]st|sachsen[- ]anhalt|^halle$|zeitz", 
            location) ~"sachsen-anhalt",
      grepl("germany[- ]th|thuring.*|^gera$", 
            location) ~"thüringen",
      
      
      # Korea
      grepl("korea", 
            location) ~ "south korea",
      grepl('^gg$|gyeonggi*', 
            location) ~ "gyeonggido",
      grepl('^gw$', 
            location) ~ "gang'weondo",
      grepl('^cn$|chungcheongnam*', 
            location) ~ "chungcheongnamdo",
      grepl('^cb$|chungcheongbuk*',
            location) ~ "chungcheongbukdo",
      grepl('^gn$',
            location) ~ "gyeongsangnamdo",
      grepl('^gb$',
            location) ~ "gyeongsangbukdo",
      grepl('^jn$',
            location) ~ "jeonranamdo",
      grepl('^jb$', 
            location) ~ "jeonrabukdo",
      grepl('^jj$',
            location) ~ "jejudo",
      grepl('^sw$', 
            location) ~ "seoul teugbyeolsi", 
      grepl('^ic$',
            location) ~ "incheon gwang'yeogsi", 
      grepl('^dj$',
            location) ~ "daejeon gwang'yeogsi", 
      grepl('^sj$|sejong teugbyeolsi', 
            location) ~ "sejong teugbyeolsi", 
      grepl('^gj$',
            location) ~ "gwangju gwang'yeogsi", 
      grepl('^dg$',
            location) ~ "daegu gwang'yeogsi", 
      grepl('^ps$', 
            location) ~ "busan gwang'yeogsi", 
      grepl('^us$', 
            location) ~ "ulsan gwang'yeogsi", 
      
      
      # Netherlands
      grepl('drenthe|westerveld|borger-odoorn|tynaarlo|eext|koekange|middendrexhe|noordenveld',
            location)  ~ 'drenthe',
      grepl('flevoland|almere|lelystad|noordoostpolder|noordosterpolder|marker wadden|stroe waddenkust',
            location) ~ 'flevoland',
      grepl('zeewolde tulpeiland|zuidlaarmeergebiede oost polder',
            location) ~ 'flevoland',
      grepl('friesland|fryslan|tytsjerksteradiel|eastermar|heerenveen|schiermonnikoog|^ameland$', 
            location) ~ "friesland",
      grepl('kollum|^vlieland.*|ameland|noord[[:space:]]{0,1}buren waddenkust', 
            location) ~ "friesland",
      grepl('gelderland|zutphen|zevenaar|westvoort|westervoort|andelst|arnhem|bennekom|deliemers',
            location) ~ 'gelderland',
      grepl('spijk|doetichem|doetinchem|epe|ermelo|gelmonde|ingen|klarenbeek|lochem|oldenbroek|huissen',
            location) ~ 'gelderland',
      grepl('putten|zwolle bergkloosterweg|apeldoorn boogschutterstraat',
            location) ~ 'gelderland',
      grepl('groningen|zuidhorn|oldambt|spijk|lauwersmeer|haren boeremapark', 
            location)  ~ 'groningen',
      grepl('limburg|venlo|gennep|kerkrade|landgraaf|ottersum|maastricht',
            location) ~ 'limburg',
      grepl('noord brabant|wernhout|vlijmen|beekendonk|best|boekel|uden|agatha|schaik|rosmalen|^reek&',
            location) ~ 'noord-brabant',
      grepl('lierop|oss|north brabant|oirschot|lithse eendenkooi',
            location) ~ 'noord-brabant',
      grepl('noord holland|bloemendaal|heerhugowaard|hilversum|huizen|naarden|reek|texel de schorren', 
            location) ~ 'noord-holland',
      grepl('wie[r]{0,1}ingerwerf|den oever|hippolytushoef|kreupel|onderdijk|^obdam$', 
            location) ~ 'noord-holland',
      grepl('ijmuiden forteiland|zwaagdijk.*|schellinkhout zuiderdijk|schagen schagerweg', 
            location) ~ 'noord-holland',
      grepl('overijssel|wierden|enschede|hardenberg|heino|losser|overdinkel|raalte',
            location)  ~ 'overijssel',
      grepl('utrecht|bilthoven|bosch en duin|bunnik|debilt|derondevenen|soest|haarzuilens|houten', 
            location)  ~ 'utrecht',
      grepl('ijsselstein|langbroek|nieuwegein|ijsselstein|woerden|leusden|vleuten.*|portengen|eemmeer dode hond', 
            location)  ~ 'utrecht',
      grepl('zeeland|terneuzen|reimerswaal|grenspad|middelburg|weversinlaag|weversinlaag|neeltje jans', 
            location) ~ 'zeeland',
      grepl('walcheren t vroon westkapelle|sloehaven', 
            location) ~ 'zeeland',
      grepl('zuid holland|zoeterwoude|westland|den haag|rotterdam|rijswijk|leiden|oud albas|oud alblas',
            location) ~ 'zuid-holland',
      grepl('nederlek|hardinxveldgiessendam|south holland|eendenkooi stolwijk|haringvliet|bliek$', 
            location) ~ 'zuid-holland',
      grepl('kwade hoek stellendam|stellendam scheelhoekeiland', 
            location) ~ 'zuid-holland',
      
      
      # Nigeria
      grepl("abuja",
            location) ~"abuja capital territory",
      grepl("abia state", 
            location) ~"abia",
      grepl("adamawa state",
            location) ~"adamawa",
      grepl("akwa ibom state", 
            location) ~"akwa ibom",
      grepl("aguata lga|anambra state",
            location) ~"anambra",
      grepl("tilden fulani toro lga|bauchi state", 
            location) ~"bauchi",
      grepl("bayelsa state", 
            location) ~"bayelsa",
      grepl("benue state", 
            location) ~"benue",
      grepl("borno state", 
            location) ~"borno",
      grepl("corss river state", 
            location) ~"cross river",
      grepl("delta state",
            location) ~"delta",
      grepl("ebonyi state", 
            location) ~"ebonyi",
      grepl("benin edo|edo state", 
            location) ~"edo",
      grepl("ekiti state",
            location) ~"ekiti",
      grepl("enugu state", 
            location) ~"enugu",
      grepl("gombe state|pantaini lbm|^akko[ gombe]{0,1}", 
            location) ~"gombe",
      grepl("imo state", 
            location) ~"imo",
      grepl("dutse lga|jigawa state",
            location) ~"jigawa",
      grepl("goni gora chikun|kaduna state", 
            location) ~"kaduna",
      grepl("kano state|gwale lga|gunduwawa gezawa|dawakin tofa",
            location) ~"kano",
      grepl("katsina state", 
            location) ~"katsina",
      grepl("kebbi state",
            location) ~"kebbi",
      grepl("kogi state", 
            location) ~"kogi",
      grepl("kwara state",
            location) ~"kwara",
      grepl("lagos state", 
            location) ~"lagos",
      grepl("keffi lgc|nasarawa state",
            location) ~"nassarawa",
      grepl("niger state", 
            location) ~"niger",
      grepl("ogun state",
            location) ~"ogun",
      grepl("ondo state", 
            location) ~"ondo",
      grepl("osun state", 
            location) ~"osun",
      grepl("oyo state", 
            location) ~"oyo",
      grepl("plateau state|jos north|maraban jos|rayfield jos south", 
            location) ~"plateau",
      grepl("port harcourt|rivers state|nkpor obio nkpor lgc", 
            location) ~"rivers",
      grepl("sokoto state", 
            location) ~"sokoto",
      grepl("taraba state",
            location) ~"taraba",
      grepl("yobe state",
            location) ~"yobe",
      grepl("zamfara state", 
            location) ~"zamfara",
      
      grepl('bosnia and herzegovina', location) ~ 'bosnia-herzegovina',
      
      grepl('sva |^sva$|north america|zeebrugge belgi', location) ~ 'NA',
      
      # Indonesia
      grepl("hulu sungai utara|banjarbaru", location) ~ "kalimantan selatan",
      grepl("republic of georgia", location) ~ "georgia",
      grepl("czech republic|czechia", location) ~ "czech republic",
      # grepl('laos', location) ~ "lao people's democratic republic",
      
      .default = location
    )) %>%
    
    rowwise() %>%
    mutate(loc = getLocation(location))
    as_tibble() %>%
    unnest(loc)  %>%
    
    
    #########################################################################################
  # Format source before taxa allocation
  mutate(source = gsub("_", " ", tolower(source))) %>%
    mutate(source = str_trim(source)) %>%
    mutate(source= gsub('gray', 'grey', source)) %>%
    mutate(primary_com_name = case_when(
      
      grepl("spot-billed duck", source) ~ "eastern spot-billed duck",
      grepl("crested grebe", source) ~ "great crested grebe",
      grepl("eastern curlew", source) ~ "far eastern curlew",
      grepl("bean goose", source) ~ "taiga/tundra bean-goose",
      grepl("surf scooter", source) ~ "surf scoter",
      grepl("northern goshawk|^goshawk$", source) ~ "eurasian goshawk", # add & for country
      grepl("coopers s hawk", source) ~ "cooper's hawk",
      grepl("rosss goose", source) ~ "ross's goose",
      grepl("madarin duck", source) ~ "mandarin duck",
      grepl("bar headed goose", source) ~ "bar-headed goose",
      grepl('grey goose|graylag goose', source) ~ 'greylag goose',
      grepl("shoveler", source) ~ "northern shoveler",
      grepl("white-fronted goose", source) ~ "greater/lesser white-fronted goose",
      grepl("mallard duck", source) ~ "mallard",
      grepl("chicken|layer hen|laying hen|broiler|poultry|chichen|^hen$|^rooster$|^layer$", source) ~ "red junglefowl (domestic type)",
      grepl("domestic goose|pomeranian goose|embden goose|sebastopol goose|anser anser domesticus|american buff goose|african goose|rural goose", source) ~ "domestic goose sp. (domestic type)",
      grepl("domestic duck|pekin duck|mule duck|runner duck|cascade duck|rural duck", source) ~ "mallard (domestic type)",
      grepl('^buzzard$', source) ~ 'common buzzard',
      grepl('knot wader', source) ~ 'red knot',
      grepl('brent goose|brant goose', source) ~ 'brant',
      grepl('eagle owl|eurasian[- ]{0,1}eagle[- ]{0,1}owl', source) ~'eurasian eagle-owl',
      grepl('canade goose|canada goose', source) ~ 'canada goose',
      grepl('european herring gull', source) ~ 'herring gull',
      grepl('towny owel', source) ~ 'tawny owl',
      grepl('gadwall duck', source) ~ 'gadwall',
      grepl('lesser snow goose', source) ~ 'snow goose',
      grepl("pink footed goose", source) ~ "pink-footed goose",
      grepl("european wigeon", source) ~ "eurasian wigeon",
      grepl('western jackdaw', source) ~ 'eurasian jackdaw',
      grepl('spotbill duck', source) ~ "indian spot-billed duck",
      grepl("whistling duck", source) ~ "whistling-duck sp.",
      grepl('pink footed goose', source) ~ 'pink-footed goose',
      grepl('falcated teal', source) ~ 'falcated duck',
      grepl('^pintail$', source) ~ 'northern pintail',
      grepl('grey plover', source)~ "black-bellied plover",
      grepl('greenwing duck', source) ~ "green-winged teal",
      grepl('franklins gull', source) ~ "franklin's gull",
      grepl("hartlaubs gull", source) ~ "hartlaub's gull",
      grepl('long-tailed skua', source) ~ "long-tailed jaeger",
      grepl('amazon parrot', source) ~ "amazona sp.",
      grepl('catalina macaw', source) ~ "large macaw sp.",
      grepl('coopers{0,2} hawk', source) ~ "cooper's hawk",
      grepl("bufflehead", source) ~ "bufflehead",
      grepl("great-white pelican", source) ~ "great white pelican",
      grepl("nene goose", source) ~ "hawaiian goose",
      grepl("american pelican", source) ~ "american white pelican",
      grepl("white breasted cormorant", source) ~ "great cormorant",
      grepl("northern giant petrel", source) ~ "northern giant-petrel",
      grepl("african fish eagle", source) ~ "african fish-eagle",
      grepl("european white stork", source) ~ "white stork",
      grepl("nothern gannet", source) ~ "northern gannet",
      grepl("chinese ringneck pheasant", source) ~ "ring-necked pheasant",
      grepl("carolina duck", source) ~ "wood duck",
      grepl("african black oystercatcher", source) ~ "african oystercatcher",
      grepl("red-backed-hawk", source) ~ "variable hawk",
      grepl("black-headead gull", source) ~ "black-headed gull",
      grepl("western screech owl", source) ~ "western screech-owl", 
      grepl("harris{1,2} hawk", source) ~ "harris's hawk",
      grepl("northwestern crow", source) ~ "american crow",    
      grepl( "lady amhersts pheasant", source) ~ "lady amherst's pheasant",
      grepl("barn owl", source) ~ 'barn owl',
      grepl("arc{0,1}tic tern", source) ~ "arctic tern",
      grepl("black-swan", source) ~ "black swan",
      grepl("swainsons hawk", source ) ~"swainson's hawk",
      grepl('^reed warbler$', source) ~ "acrocephalus sp.",
      grepl("common|eurasian teal|green[- ]{0,1}winged[- ]{0,1}teal", source) ~ "green-winged teal", 
      grepl("^oystercatcher$", source) ~ "oystercatcher sp.", 
      grepl("^curlew$", source) ~ "curlew sp.", 
      grepl("blue-winged teal", source) ~ "blue-winged teal", 
      grepl("swift tern", source) ~ "great crested tern",   
      grepl("shelduck" , source) ~ "dabbling duck sp.",
      grepl("rt-hawk", source) ~  "red-tailed hawk",
      grepl("sea eagle", source) ~ "eagle sp.",
      grepl("^moorhen$", source) ~ "rail/crake sp.",
      grepl("^kittiwake$", source) ~ "black-legged/red-legged kittiwake",
      grepl('mute[- ]{0,1}swan', source) ~"mute swan",
      grepl('^rhea$', source)  ~ "rhea sp.",
      
      
      # latin names
      grepl("anas crecca", source) ~ "green-winged teal",
      grepl('tyto alba', source) ~ 'barn owl',
      grepl('gallus gallus', source) ~ 'red junglefowl (domestic type)',
      grepl("anas platyrhynchos", source) ~ "mallard",
      grepl("anser albifrons", source) ~ "greater white-fronted goose",
      grepl("anser brachyrhynchus", source) ~ "pink-footed goose",
      grepl("anser anser", source) ~ "greylag goose",
      grepl("anser fabalis", source) ~ "taiga/tundra bean-goose",
      grepl("chlidonias hybrida", source) ~ "whiskered tern",
      grepl("cygnus columbianus", source) ~ "tundra swan",
      grepl("falco peregrinus", source) ~ "peregrine falcon",
      grepl("anser fabalis", source) ~ "taiga/tundra bean-goose",
      grepl('larus argentatus', source) ~ 'herring gull',
      grepl('pica pica', source) ~ 'eurasian magpie',
      grepl("ciconia ciconia", source) ~ 'white stork',
      grepl("buteo buteo", source) ~ 'common buzzard',
      grepl('podiceps cristatus', source) ~ 'great crested grebe',
      grepl('anas poecilorhyncha|anas poecilohyncha', source) ~ "indian spot-billed duck",
      grepl('garrulus glandarius', source) ~ 'eurasian jay',
      grepl('columba palumbus', source) ~ "common wood-pigeon",
      grepl('branta canadensis', source) ~ 'canada goose',
      grepl('cygnus olor', source) ~"mute swan",
      grepl('spatula clypeata', source) ~ "northern shoveler",
      grepl('tadorna tadorna', source) ~ "common shelduck",
      grepl('anseriformes', source) ~ 'waterfowl sp.',
      grepl('sterna hirundo', source) ~ "common tern",
      grepl("meleagris gallopavo", source) ~ "wild turkey",
      grepl("branta leucopsis", source) ~ "barnacle goose",
      grepl("alopochen aegyptiaca", source) ~ "egyptian goose",
      grepl("numenius arquata", source) ~ 'eurasian curlew',
      grepl("pavo cristatus", source) ~ "indian peafowl",
      grepl("pelecanus", source) ~ "pelican sp.",             
      grepl("corvus corvus", source) ~ "crow sp.",              
      grepl("cygnus cygnus", source) ~ "whooper swan",
      grepl("chroicocephalus ridibundus", source) ~ "black-headed gull",
      grepl("otus scops", source) ~ "eurasian scops-owl",                 
      grepl("aythya fuligula", source) ~ "tufted duck",             
      grepl("corvus monedula", source) ~ "eurasian jackdaw",            
      grepl("morus bassanus", source)~ "northern gannet", 
      grepl("tachybaptus ruficollis", source) ~ "little grebe",    
      grepl("egretta garzetta", source)~"little egret",            
      grepl("pelecanus conspicillatus", source) ~ "australian pelican",   
      grepl("phasianus colchicus", source) ~ "ring-necked pheasant",   
      grepl("gypaetus barbatus", source) ~ "bearded vulture",
      grepl("pelecanus occidentalis", source) ~ "brown pelican",
      grepl("ardea cinerea", source) ~ "grey heron",
      grepl("thalasseus sandvicensis", source) ~ "sandwich tern",    
      grepl("phalacrocorax carbo", source) ~ "great cormorant",
      grepl("laridae", source) ~ "gull/tern sp.",
      grepl("phasianidae", source) ~ "pheasant sp.",
      grepl('anatidae', source) ~ "teal sp.",    
      grepl('bubo bubo', source) ~'eurasian eagle-owl',
      grepl('larus canus', source) ~ 'common gull',
      
      
      # generics
      grepl("^swan$", source) ~ "swan sp.",
      grepl("^cormorant$", source) ~ 'cormorant sp.',
      grepl("^egret$", source) ~ "white egret sp.",
      grepl("^crane$|black-throated crane", source) ~ "crane sp.",
      grepl("^goose$|wild geese|wild goose", source) ~ "goose sp.",
      grepl('^stork$', source) ~ 'stork sp.',
      grepl('^ostrich$' ,source) ~ "common ostrich",
      grepl("^egret$", source) ~ "egret sp.",
      grepl('lapwing', source) ~ 'lapwing sp.',
      grepl('^wigeon$', source) ~ "eurasian/american wigeon",
      grepl('^pochard$', source) ~ "aythya sp.",
      grepl('coot', source) ~ 'coot sp.',
      grepl("waterfowl", source) ~ "waterfowl sp.",
      grepl("peacock|peafowl", source) ~ "indian peafowl (domestic type)",
      
      grepl("shorebird|sandpiper", source) ~ "shorebird sp.",
      grepl("black-backed gull|^gull$|seagull|seabird", source) ~ "gull sp.",
      grepl("^crow$|^raven$|american raven", source) ~ "crow/raven sp.",
      grepl("^duck$|wild duck|migratory duck", source) ~ "duck sp.",
      grepl("^grebe$", source) ~ "grebe sp.",
      grepl("^owl$", source) ~ "owl sp.",
      grepl('^hawk$', source) ~ "hawk sp.",
      grepl('^eagle$', source) ~ 'eagle sp.',
      grepl('steamer duck', source) ~ "steamer-duck sp.",
      grepl('^falcon$', source) ~ "falcon sp.",
      grepl("^flamingo$", source) ~ "flamingo sp.",
      grepl("^pelican$", source) ~ "pelican sp.",
      grepl("^pheasant$|^partridge$", source) ~ "pheasant sp.",
      grepl("^pigeon$", source) ~ "pigeon/dove sp.",
      grepl("^teal$", source) ~ "teal sp.",
      grepl('^ibis$', source) ~ 'ibis sp.',
      grepl('^rhea$', source) ~ 'rhea sp.',
      
      
      grepl('^gannet$', source) ~ "sulid sp.",
      grepl("^heron$", source) ~"heron sp.",
      grepl("^aquatic bird$|^avian$|^bird$|wild[- ]{0,1}birds{0,1}|backyard bird", source) ~ "bird sp.",
      source %in% c("en", "env", "environment", "enviroment", "environmental", 'environment sample',
                    "water", "wild bird feces") ~ "environment",
      
      # Mammals
      grepl("red fox|^fox$|vulpes vulpes", source) ~ "red fox",
      grepl("r{0,1}attus norvegicus", source) ~ "brown rat",
      grepl('^seal$|sea lion', source) ~ 'seal sp.', 
      grepl('^mink$|wild mink', source) ~ 'mink sp.',
      grepl("harbou{0,1}r seal", source) ~ 'harbour seal',
      grepl("mustela furo", source) ~"ferret",
      grepl("mustela putorius", source) ~"european polecat",
      grepl("stone marten", source) ~"beech marten",
      grepl("bear$", source) ~"bear sp.",
      grepl("dolphin$|porpoise", source) ~"cetacea sp.",
      grepl("pekania pennanti", source) ~"fisher",
      grepl("^fox$", source) ~"fox sp.",
      grepl("^lynx$", source) ~"lynx sp.",
      grepl("^polecat$", source) ~"mustelid sp.",
      grepl('^skunk$', source) ~ 'skunk sp.',
      grepl('^otter$', source) ~ 'otter sp.',
      grepl('^lion$|^cat$|domestic cat|feline', source) ~ 'feline sp.',
      grepl('canine', source) ~ 'canine sp.',
      
      # locations
      grepl('sparrowhawk', source) & grepl('europe', region) ~ 'eurasian sparrowhawk',
      grepl("^magpie$", source) & grepl('europe', region) ~ "eurasian magpie",
      grepl("^magpie$", source) & grepl('eastern asia', region) ~ "oriental/eurasian magpie",
      grepl("^magpie$", source) & grepl('northern america', region) ~ "black-billed magpie",
      grepl("^kestrel$", source) & grepl('northern america', region) ~ "american kestrel",
      grepl("^kestrel$", source) & grepl('europe', region) ~ "lesser/eurasian kestrel",
      grepl("^guinea {0,1}fowl", source) & grepl('africa', region) ~ "crested guineafowl sp.",
      grepl('^guinea {0,1}fowl', source) & !grepl('africa', region) ~ "helmeted guineafowl (domestic type)",
      grepl('^fulmar$', source) & !grepl('south america', region) ~ "northern fulmar",
      grepl('^fulmar$', source) & grepl('south america', region) ~ "southern fulmar",
      grepl('^gannet$', source) & grepl('europe', region) ~ "northern gannet",
      grepl('jungle crow', source) ~ "crow/raven sp.",
      grepl('oystercatcher', source) & grepl('europe', region) ~ 'eurasian oystercatcher',
      grepl("turkey|^pavo$", source) & grepl('central america', region) ~ "wild turkey", 
      grepl("turkey|^pavo$", source) & !grepl('central america', region) ~ "wild turkey (domestic type)", 
      grepl("turnstone", source) & grepl('europe', region) ~ "ruddy turnstone",
      source == 'sea eagle' & grepl('europe', region) ~"white-tailed eagle",
      grepl("quail", source) & !grepl('america', region) ~ "old world quail sp.", 
      grepl("quail", source) & grepl('america', region) ~ "new world quail sp.", 
      grepl('^vulture$', source) &  grepl('america', region) ~ "new world vulture sp.",
      grepl('^vulture$', source) & !grepl('america', region) ~ "old world vulture sp.",
      grepl('^otter$', source) & grepl('europe', region) ~ 'eurasian otter',
      grepl('badger', source) & grepl('europe', region) ~ 'european badger',
      
      .default = source 
    )) %>%
    mutate(is.problem.bird = case_when(primary_com_name %in% tolower(birds$PRIMARY_COM_NAME) ~ FALSE, 
                                       .default = TRUE))%>%
    
    # Get host taxonomy
    rowwise() %>%
    mutate(taxonomy = getTaxonomyForName(primary_com_name))%>%
    as_tibble() %>%
    unnest(taxonomy) %>%
    rename_with(.fn = ~ tolower(.x)) %>%
    
    # binary host information
    mutate(is.domestic = case_when(grepl("domestic type", sci.name) ~ "domestic",
                                   .default = "wild"
    )) %>%
    mutate(is.bird = case_when(grepl("ormes$", order) ~ "bird",
                               .default = "other"
    )) %>%
    
    # Format column names
    #dplyr::select(-c(virus_species, id_unsure)) %>%
    dplyr::rename(
      virus.subtype = subtype,
      collection.date = date,
      collection.datedecimal = decimal.date,
      collection.datemonth = month.date,
      collection.region.name = region,
      collection.country.name = name_country,
      collection.country.code = code_country,
      collection.country.lat = lat_country,
      collection.country.long = long_country,
      collection.subdiv1.name = name_subdiv1,
      collection.subdiv1.code = code_subdiv1,
      collection.subdiv1.lat = lat_subdiv1,
      collection.subdiv1.long = long_subdiv1,
      collection.subdiv2.name = name_subdiv2,
      collection.subdiv2.code = code_subdiv2,
      collection.subdiv2.lat = lat_subdiv2,
      collection.subdiv2.long = long_subdiv2,
      host.order = order,
      host.family = family,
      host.sciname = sci.name,
      host.commonname = primary_com_name,
      host.isbird = is.bird,
      host.isdomestic = is.domestic
    ) %>%
    dplyr::relocate(
      virus.subtype,
      isolate.id,
      isolate.name,
      collection.date,
      collection.datedecimal,
      collection.datemonth,
      host.order,
      host.family,
      host.sciname,
      host.commonname,
      host.isbird,
      host.isdomestic,
      collection.region.name,
      collection.country.name,
      collection.country.code,
      collection.country.lat,
      collection.country.long,
      collection.subdiv1.name,
      collection.subdiv1.code,
      collection.subdiv1.lat,
      collection.subdiv1.long,
      collection.subdiv2.name,
      collection.subdiv2.code,
      collection.subdiv2.lat,
      collection.subdiv2.long
    ) %>%
      
    # Select columns
    select(-c(source, location)) %>%
    mutate(across(everything(), .fns = ~ gsub('^NA$', NA, .x)))
  
  return(data)
  
}


MergeReassortantData <- function(data, newdata){
  stopifnot(any('isolate_id' %in% colnames(newdata)))
  
  joineddata <- data %>% 
    left_join(newdata, 
              by = join_by(isolate.id)) %>%
    
    # Resolve discrepancies
    mutate(collection.datedecimal = case_when(
      collection.datedecimal.x != collection.datedecimal.y ~ collection.datedecimal.y
    ))
    
    
    
  
  
    newdata %>% 
    select(c(isolate_id, location, contains('cluster'), date_year_month)) %>%
    separate_wider_delim(location, 
                         delim = ' / ', 
                         names = c('region',
                                   'collection.country.name', 
                                   'collection.subdiv1.name',
                                   'collection.subdiv2.name',
                                   'collection.subdiv3.name'),
                         too_few = 'align_start') %>%
    rename_with(.fn = ~gsub('_', '.', .x))
  
  # Get reassortant data from metadata
  out <- left_join(data, 
                   newdata,
                   by = join_by(isolate.id)) %>%
    
    # Get missing dates from metadata csv
    mutate(collection.dateweek = dplyr::coalesce(collection.dateweek, 
                                                 date.year.month)) %>%
    
    mutate(collection.tipdate = case_when(is.na(collection.date) ~ collection.dateweek,
                                          .default = as.character(collection.date))) #%>%
  # select(-c(date.month.year, month.date))

  return()

}
