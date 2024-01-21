FormatMetadata <- function(data){
  data <- data %>%
    # Format source before taxa allocation
    mutate(source = gsub("_", " ", tolower(source))) %>%
    mutate(source = str_trim(source)) %>%
    mutate(source= gsub('gray', 'grey', source)) %>%
    mutate(primary_com_name = case_when(
      grepl("common|eurasian teal|green-winged-teal", source) ~ "green-winged teal", # This one is a pain
      grepl("spot-billed duck", source) ~ "eastern spot-billed duck",
      grepl("crested grebe", source) ~ "great crested grebe",
      grepl("eastern curlew", source) ~ "far eastern curlew",
      grepl("bean goose", source) ~ "taiga/tundra bean-goose",
      grepl("surf scooter", source) ~ "surf scoter",
      grepl("northern goshawk|^goshawk$", source) ~ "eurasian goshawk", # add & for country
      grepl('sparrowhawk', source) & grepl('England', location) ~ 'eurasian sparrowhawk',
      grepl("coopers s hawk", source) ~ "cooper's hawk",
      grepl("rosss goose", source) ~ "ross's goose",
      grepl("turnstone", source) & location == "Netherlands" ~ "ruddy turnstone",
      grepl("madarin duck", source) ~ "mandarin duck",
      grepl("bar headed goose", source) ~ "bar-headed goose",
      grepl('grey goose|graylag goose', source) ~ 'greylag goose',
      grepl("shoveler", source) ~ "northern shoveler",
      grepl("white-fronted goose", source) ~ "greater/lesser white-fronted goose",
      grepl("mallard duck", source) ~ "mallard",
      grepl("chicken|layer hen|laying hen|broiler|poultry|chichen", source) ~ "red junglefowl (domestic type)",
      grepl("domestic goose|pomeranian goose", source) ~ "domestic goose sp. (domestic type)",
      grepl("domestic duck|pekin duck|mule duck", source) ~ "mallard (domestic type)",
      grepl('^buzzard$', source) ~ 'common buzzard',
      grepl('knot wader', source) ~ 'red knot',
      grepl("turkey", source) & location != "Mexico" ~ "wild turkey",
      grepl('^kestrel$', source) & grepl('Germany|England|Italy', location) ~ 'eurasian kestrel',
      grepl('brent goose', source) ~ 'brant',
      grepl('eagle owl', source) ~'eurasian eagle-owl',
      grepl('jungle crow', source) & location == 'Japan' ~ 'large-billed crow',
      grepl('^magpie$' ,source) & location == 'Idaho' ~ 'black-billed magpie',
      grepl('oystercatcher', source) & grepl('Germany', location) ~ 'eurasian oystercatcher',
      grepl('canade goose', source) ~ 'canada goose',
      grepl('european herring gull', source) ~ 'herring gull',
      grepl('towny owel', source) ~ 'tawny owl',
      grepl('gadwall duck', source) ~ 'gadwall',
      grepl('lesser snow goose', source) ~ 'snow goose',
      grepl("pink footed goose", source) ~ "pink-footed goose",
      grepl("european wigeon", source) ~ "eurasian wigeon",
      source == 'sea eagle' & grepl('Norway', location)~"white-tailed eagle",
      grepl('guinea fowl', source) & grepl('Italy|Estonia|Germany', location) ~ "helmeted guineafowl (domestic type)",
      grepl('western jackdaw', source) ~ 'eurasian jackdaw',
      grepl('^gannet$', source) & grepl('Scotland', location) ~ "northern gannet",
      grepl('spotbill duck', source) ~ "indian spot-billed duck",
      grepl("whistling duck", source) ~ "whistling-duck sp.",
      grepl('pink footed goose', source) ~ 'pink-footed goose',
      grepl('falcated teal', source) ~ 'falcated duck',
      grepl('^pintail$', source) ~ 'northern pintail',
      grepl('grey plover', source)~ "black-bellied plover",
      grepl('greenwing duck', source) ~ "green-winged teal",
      grepl('franklins gull', source) ~ "franklin's gull",
      
      
      
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
      
      grepl('^seal$', source) ~ 'seal sp.', 
      grepl('^mink$', source) ~ 'mink sp.',
      grepl('^otter$', source) & grepl('Scotland', location) ~ 'eurasian otter',
      grepl('badger', source) ~ 'european badger',
      
      
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
      grepl('^pochard$', source) ~ "aythya sp.",
      grepl('coot', source) ~ 'coot sp.',
      grepl("waterfowl", source) ~ "waterfowl sp.",
      grepl("peacock|peafowl", source) ~ "indian peafowl (domestic type)",
      grepl("quail", source) & location != "America" ~ "old world quail sp.", # requires better location
      grepl("shorebird|sandpiper", source) ~ "shorebird sp.",
      grepl("black-backed gull|^gull$|seagull", source) ~ "gull sp.",
      grepl("^crow$", source) ~ "crow/raven sp.",
      grepl("^duck$|wild duck|migratory duck", source) ~ "duck sp.",
      grepl("^grebe$", source) ~ "grebe sp.",
      grepl("^owl$", source) ~ "owl sp.",
      grepl('^hawk$', source) ~ "hawk sp.",
      grepl('^eagle$', source) ~ 'eagle sp.',
      grepl('steamer duck', source) ~ "steamer-duck sp.",
      grepl('^falcon$', source) ~ "falcon sp.",
      grepl("^pelican$", source) ~ "pelican sp.",
      grepl("^pheasant$", source) ~ "pheasant sp.",
      grepl("^pigeon$", source) ~ "pigeon/dove sp.",
      grepl("^teal$", source) ~ "teal sp.",
      grepl('^ibis$', source) ~ 'ibis sp.',
      source %in% c("aquatic bird", "avian", "bird", "wild bird", "wild bird", "backyard bird") ~ "bird sp.",
      source %in% c("en", "env", "environment", "enviroment", "environmental", 'environment sample',
                    "water", "wild bird feces") ~ "environment",
      grepl("red fox|^fox$", source) ~ "red fox",
      grepl("attus norvegicus", source) ~ "brown rat",
      .default = source
    )) %>%
    
    # Get host taxonomy
    rowwise() %>%
    mutate(taxonomy = getTaxonomyForName(primary_com_name)) %>%
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
    
    # Location (iso name and subdivision)
    #mutate(location = tolower(location)) %>%
    mutate(location = gsub("_", " ", tolower(location))) %>%
    
    mutate(location = case_when(
      
      # Korea
      grepl("korea", location) ~ "south korea",
      
      # PRC
      grepl("yunnan.*", location) ~ "yunnan",
      grepl('rongcheng', location) ~ 'shandong',
      grepl("inner mongolia|tumuji", location) ~ "nei mongol",
      grepl("hong kong|hongkong", location) ~ "xianggang (hong-kong)",
      grepl("xizang$", location) ~ "xizang",
      grepl("ningxia$", location) ~ "ningxia",
      grepl("xinjiang$", location) ~ "xinjiang",
      grepl("guangxi$|guilin|hechi", location) ~ "guangxi",
      grepl("tibet", location) ~ "xizang",
      grepl("hangzhou|nanji$", location) ~ "zhejiang",  # confirm nanji
      grepl("sanmenxia|dongting", location) ~ "henan",
      grepl("heinan", location) ~ "hainan", # confirm 
      grepl("wuhan", location) ~ "hebei",
      grepl('^fujain$', location) ~ 'fujian',
      grepl('suzhou|xuzhou', location) ~ 'jiangsu',
      grepl("changsha", location) ~ "hunan",
      grepl('tianjing', location) ~ 'tianjin',
      grepl('qingyuan|foshan', location) ~ 'guangdong',
      grepl("northern china", location) ~ "beijing", # Interpretation from paper - check with Lu
      grepl("southwestern china", location) ~ "shanxi", # Interpretation from paper - check with Lu
      grepl('eastern china', location) ~ 'zhejiang', #Inferred - if we need to be 100% sure then change to country level
      
      
      # Canada
      grepl("^ab$", location) ~ "alberta",
      grepl("^bc$", location) ~ "british columbia",
      grepl('^mb$', location) ~'manitoba',
      
      # Russian Federation
      grepl('russian federation', location) ~ 'russia',
      grepl('khabarovsk', location) ~ "khabarovskiy kray",
      grepl('chelyabinsk', location) ~ "chelyabinskaya oblast'",
      grepl("astrakhan", location) ~ "astrakhanskaya oblast'",
      grepl("novosibirsk region|chany lake|^chany$", location) ~ "novosibirskaya oblast'",
      grepl("omsk.*", location) ~ "omskaya oblast'",
      grepl("rostov-on-don", location) ~ "rostovskaya oblast'",
      grepl("russia primorje", location) ~ "primorskiy kray",
      grepl("sakhalin", location) ~ "sakhalinskaya oblast'",
      grepl("yakutia", location) ~ "sakha, respublika [yakutiya]",
      grepl("kostroma", location) ~ "kostromskaya oblast'",
      grepl("amur region", location) ~ "amurskaya oblast'",
      grepl("buryatia", location) ~ "buryatiya, respublika",
      grepl('north ossetia-alania', location) ~ 'severnaya osetiya-alaniya, respublika',
      grepl('dagestan', location) ~ 'dagestan, respublika',
      grepl('stavropol', location) ~ "stavropol'skiy kray",
      grepl('krasnodar', location) ~ "krasnodarskiy kray", 
      grepl('tyumen', location) ~ "tyumenskaya oblast'",
      grepl("magadan", location) ~ "magadanskaya oblast'",
      grepl('saratov', location) ~ "saratovskaya oblast'",
      grepl('kurgan*', location) ~ "kurganskaya oblast'",
      grepl('central russia' ,location) ~ 'russia',
      
      # USA
      grepl("north dakota", location) ~ "north dakota",
      
      # UK
      grepl("england", location) ~ "england and wales",
      # Japan
      grepl("tsukuba", location) ~ "ibaraki",
      
      # Spain
      grepl("castillalamancha", location) ~ "castilla-la mancha",
      
      # Germany 
      location == 'germany-mv' ~ 'mecklenburg-vorpommern',
      location == 'germany-nw' ~ 'nordrhein-westfalen',
      location == 'germany-be' ~ 'berlin',
      location == 'germany-by' ~ 'bayern',
      location =='germany-st' ~ 'sachsen-anhalt',
      location == 'germany-bb' ~ 'brandenburg',
      location == 'germany-ni' ~ 'niedersachsen',
      location == 'germany-sh' ~ 'schleswig-holstein',
      location == 'germany-sn' ~ 'sachsen',
      location == 'germany-he' ~ 'hessen',
      location == 'germany-hh' ~ 'hamburg',
      location == 'germany-bw' ~ 'baden-württemberg',
      location == 'germany-th' ~ 'thüringen',
      location == 'germany-rp' ~ 'rheinland-pfalz',
      
      # Korea
      grepl('^gg$|gyeonggi*', collection.subdiv1.name) ~ "gyeonggido",
      grepl('^gw$', collection.subdiv1.name) ~ "gang'weondo",
      grepl('^cn$|chungcheongnam*', collection.subdiv1.name) ~ "chungcheongnamdo",
      grepl('^cb$|chungcheongbuk*', collection.subdiv1.name) ~ "chungcheongbukdo",
      grepl('^gn$', collection.subdiv1.name) ~ "gyeongsangnamdo",
      grepl('^gb$', collection.subdiv1.name) ~ "gyeongsangbukdo",
      grepl('^jn$', collection.subdiv1.name) ~ "jeonranamdo",
      grepl('^jb$', collection.subdiv1.name) ~ "jeonrabukdo",
      grepl('^jj$', collection.subdiv1.name) ~ "jejudo",
      grepl('^sw$', collection.subdiv1.name) ~ "seoul teugbyeolsi", 
      grepl('^ic$', collection.subdiv1.name) ~ "incheon gwang'yeogsi", 
      grepl('^dj$', collection.subdiv1.name) ~ "daejeon gwang'yeogsi", 
      grepl('^sj$', collection.subdiv1.name) ~ "sejong teugbyeolsi", 
      grepl('^gj$', collection.subdiv1.name) ~ "gwangju gwang'yeogsi", 
      grepl('^dg$', collection.subdiv1.name) ~ "daegu gwang'yeogsi", 
      grepl('^ps$', collection.subdiv1.name) ~ "busan gwang'yeogsi", 
      grepl('^us$', collection.subdiv1.name) ~ "ulsan gwang'yeogsi", 
      
     
      
      
      grepl('bosnia and herzegovina', location) ~ 'bosnia-herzegovina',
      
      grepl('sva[:0-9:]', location) ~ 'NA',
      
      # Indonesia
      grepl("hulu sungai utara|banjarbaru", location) ~ "kalimantan selatan",
      grepl("republic of georgia", location) ~ "georgia",
      grepl("czech republic|czechia", location) ~ "czech republic",
      # grepl('laos', location) ~ "lao people's democratic republic",
      grepl("vietnam", location) ~ "viet nam",
      grepl("quang ninh" , location) ~  "quảng ninh" ,
      grepl('tinh lang son', collection.subdiv1.name) ~ "lạng sơn",
      
      .default = location
    )) %>%
    rowwise() %>%
    mutate(loc = getLocation(location)) %>%
    as_tibble() %>%
    unnest(loc) %>%
    
    # Format column names
    dplyr::select(-c(virus_species, id_unsure)) %>%
    dplyr::rename(
      virus.subtype = subtype,
      collection.date = date,
      collection.datedecimal = decimal.date,
      collection.dateweek = week.date,
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
      collection.dateweek,
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
    select(-c(source, location)) 
  
  return(data)
  
}
reassortant_metadata <- read_csv('./data/metadata/h5_metadata.csv')
test_reassortant <- reassortant_metadata %>%
  rename(joint_location = location) %>%
  mutate(location = coalesce(location_3,location_2)) %>%
  mutate(date = parsedate::parse_date(date) %>%
           as.Date()) %>%
  mutate(decimal.date = format(round(decimal_date(date), 2), 
                               nsmall = 2) ) %>%
  mutate(decimal.date = suppressWarnings(as.double(decimal.date))) %>%
  mutate(week.date = format(date, "%Y-%V")) %>%
  mutate(month.date = format(date, "%Y-%m")) %>%
  select(-c('joint_location', 'date_yday', 'date_frac', 'date_year', 'date_year_month', 'date_yday', 'date_year_week')) %>%
  rename_with(.fn = ~gsub('_', '.', .x)) %>%
  rename(source=host) %>%
  
  # Location (iso name and subdivision)
  #mutate(location = tolower(location)) %>%
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
    grepl("brussels-capital region|region de bruxelles capitale",
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
    grepl("krong phnum penh", location) ~"phnom penh",
    grepl("khett kandal", location) ~"kandal",
    grepl("khett prey veng", location) ~"prey veaeng",
    grepl("takeo", location) ~"taakaev",
    
    # Chile
    grepl("region de antofagasta", location) ~ "antofagasta",
    
    # Colombia
    grepl("departamento de bolivar", location) ~"bolívar",
    grepl("departamento del choco", location) ~"chocó",
    grepl("departamento de cordoba", location) ~"córdoba",
    grepl("departamento del magdalena", location) ~"magdalena",
    
    # Costa Rica
    grepl("provincia de limon", location) ~ "limón",
    
    # Croatia
    
    
    # Denmark
    
    # Ecuador
    
    # Honduras
    
    # Hungary
    
    # Iceland
    
    # Indonesia
    
    # Iraq
    
    # Israel
    
    # Kazakhstan
    
    # Latvia
    
    # Lithuania
    
    # Luxembourg
    
    # Mali
    
    # Montenegro
    
    # Netherlands
    
    # Niger
    
    # Norway
    
    # Peru
    
    # Portugal
    
    # Romania
    
    # Senegal
    grepl('region de thies',
          location) ~ "thiès",
    grepl('djoudj bird sanctuary',
          location) ~ "saint-louis",
    
    # Slovakia
    
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
    grepl("champagne-ardenne",
          location) ~"champagne-ardenne",
    grepl("corsica", 
          location) ~"corse",
    grepl("franche-comté",
          location) ~"franche-comté",
    grepl("haute-normandie", 
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
    
    # Kosovo
    

    
    
    # Russian Federation
    grepl('russian federation', 
          location) ~ 'russia',
    grepl('khabarovsk', 
          location) ~ "khabarovskiy kray",
    grepl('chelyabinsk', 
          location) ~ "chelyabinskaya oblast'",
    grepl("astrakhan", 
          location) ~ "astrakhanskaya oblast'",
    grepl("novosibirsk region|chany lake|^chany$", 
          location) ~ "novosibirskaya oblast'",
    grepl("omsk.*", 
          location) ~ "omskaya oblast'",
    grepl("rostov-on-don", 
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
    grepl('central russia' 
          ,location) ~ 'russia',
    
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
    grepl("england",
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
    grepl("lower silesian voivodeship", location) ~"dolnośląskie",
    grepl("kuyavian pomeranian voivodeship", location) ~"kujawsko-pomorskie",
    grepl("lublin voivodeship", location) ~"lubelskie",
    grepl("lubusz voivodeship", location) ~"lubuskie",
    grepl("lodzkie|odz voivodeship", location) ~"łódzkie",
    grepl("lesser poland voivodeship", location) ~"małopolskie",
    grepl("masovian voivodeship", location) ~"mazowieckie",
    grepl("opole voivodeship", location) ~"opolskie",
    grepl("pomeranian voivodeship", location) ~"pomorskie",
    grepl("silesian voivodeship", location) ~"śląskie",
    grepl("swietokrzyskie voivodeship", location) ~"świętokrzyskie",
    grepl("warmi sko mazurskie", location) ~"warmińsko-mazurskie",
    grepl("greater poland voivodeship", location) ~"wielkopolskie",
    grepl("west pomeranian voivodeship", location) ~"zachodniopomorskie",
    
    
    # Sweden
    grepl("blekinge lan", location) ~"blekinge län",
    grepl("dalarnas lan", location) ~"dalarnas län",
    grepl("gotlands lan", location) ~"gotlands län",
    grepl("hallands lan", location) ~"hallands län",
    grepl("jonkopings lan", location) ~"jönköpings län",
    grepl("kalmar lan", location) ~"kalmar län",
    grepl("skane lan", location) ~"skåne län",
    grepl("stockholms lan", location) ~"stockholms län",
    grepl("sodermanlands lan", location) ~"södermanlands län",
    grepl("uppsala lan", location) ~"uppsala län",
    grepl("vastra gotalands lan", location) ~"västra götalands län",
    grepl("ostergotlands lan", location) ~"östergötlands län",
    
    
    # Germany 
    grepl("germany", 
          location) ~"germany",
    grepl("germany-bw|baden wuerttemberg|freiburg im breisgau|rosengarten|stutensee",
          location) ~"baden-württemberg",
    grepl("germany-by|bavaria|munich|nuremberg|wuerzburg", 
          location) ~"bayern",
    grepl("germany-hb", 
          location) ~"bremen",
    grepl("germany-hh", 
          location) ~"hamburg",
    grepl("germany-he|hessen|freigericht|giessen|karben|ranstadt", 
          location) ~"hessen",
    grepl("germany-ni|lower saxony|brietlingen|dinklage|hannover|lueneburg|nordhorn|osnabruck",
          location) ~"niedersachsen", 
    grepl("osterholz-scharnbeck|theene|wingst",
          location) ~"niedersachsen",
    grepl("germany-nw|north rhine westphalia|aachen|brueggen", 
          location) ~"nordrhein-westfalen",
    grepl("germany-rp|rhineland palatinate", 
          location) ~"rheinland-pfalz",
    grepl("germany-sl", 
          location) ~"saarland",
    grepl("germany-sh|schleswig holstein|aumuehle|brickeln|luebeck|pinneberg|prohn", 
          location) ~"schleswig-holstein",
    grepl("germany-be",
          location) ~"berlin",
    grepl("germany-bb", 
          location) ~"brandenburg",
    grepl("germany-mv|mecklenburg vorpommern|grevesmuehlen|parchim|ruegen|warin", 
          location) ~"mecklenburg-vorpommern",
    grepl("germany-sn|sachsen|doberschuetz|dresden|gross dueben|leipzig|saxony",
          location) ~"sachsen",
    grepl("germany-st|sachsen-anhalt|^halle$|zeitz", 
          location) ~"sachsen-anhalt",
    grepl("germany-th|thuring.*|^gera$", 
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
    grepl('^sj$', 
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
    ########################################################################################
    
    #########################################################################################
    grepl('bosnia and herzegovina', location) ~ 'bosnia-herzegovina',
    
    grepl('sva[:0-9:]', location) ~ 'NA',
    
    # Indonesia
    grepl("hulu sungai utara|banjarbaru", location) ~ "kalimantan selatan",
    grepl("republic of georgia", location) ~ "georgia",
    grepl("czech republic|czechia", location) ~ "czech republic",
    # grepl('laos', location) ~ "lao people's democratic republic",

    
    .default = location
  )) %>%
  mutate(is.problem = case_when(location %in% tolower(geographicalmetadata$name_subdiv1) ~ FALSE, 
                                location %in% tolower(geographicalmetadata$name_country) ~ FALSE,
                                location %in% tolower(geographicalmetadata$name_subdiv2) ~ FALSE,
                                .default = TRUE))

probs <- test_reassortant %>%
  filter(is.problem == TRUE) %>%
  select(location, location.2) %>%
  distinct()

probs %>% filter(location.2 == 'Czech Republic') %>% pull(location) %>% unique()



  rowwise() %>%
  mutate(loc = getLocation(location)) %>%
  as_tibble() %>%
  unnest(loc)


  # Format source before taxa allocation
  mutate(source = gsub("_", " ", tolower(source))) %>%
  mutate(source = str_trim(source)) %>%
  mutate(source= gsub('gray', 'grey', source)) %>%
  mutate(primary_com_name = case_when(
    grepl("common|eurasian teal|green-winged-teal", source) ~ "green-winged teal", # This one is a pain
    grepl("spot-billed duck", source) ~ "eastern spot-billed duck",
    grepl("crested grebe", source) ~ "great crested grebe",
    grepl("eastern curlew", source) ~ "far eastern curlew",
    grepl("bean goose", source) ~ "taiga/tundra bean-goose",
    grepl("surf scooter", source) ~ "surf scoter",
    grepl("northern goshawk|^goshawk$", source) ~ "eurasian goshawk", # add & for country
    grepl('sparrowhawk', source) & grepl('England', location) ~ 'eurasian sparrowhawk',
    grepl("coopers s hawk", source) ~ "cooper's hawk",
    grepl("rosss goose", source) ~ "ross's goose",
    grepl("turnstone", source) & location == "Netherlands" ~ "ruddy turnstone",
    grepl("madarin duck", source) ~ "mandarin duck",
    grepl("bar headed goose", source) ~ "bar-headed goose",
    grepl('grey goose|graylag goose', source) ~ 'greylag goose',
    grepl("shoveler", source) ~ "northern shoveler",
    grepl("white-fronted goose", source) ~ "greater/lesser white-fronted goose",
    grepl("mallard duck", source) ~ "mallard",
    grepl("chicken|layer hen|laying hen|broiler|poultry|chichen", source) ~ "red junglefowl (domestic type)",
    grepl("domestic goose|pomeranian goose", source) ~ "domestic goose sp. (domestic type)",
    grepl("domestic duck|pekin duck|mule duck", source) ~ "mallard (domestic type)",
    grepl('^buzzard$', source) ~ 'common buzzard',
    grepl('knot wader', source) ~ 'red knot',
    grepl("turkey", source) & location != "Mexico" ~ "wild turkey",
    grepl('^kestrel$', source) & grepl('Germany|England|Italy', location) ~ 'eurasian kestrel',
    grepl('brent goose', source) ~ 'brant',
    grepl('eagle owl', source) ~'eurasian eagle-owl',
    grepl('jungle crow', source) & location == 'Japan' ~ 'large-billed crow',
    grepl('^magpie$' ,source) & location == 'Idaho' ~ 'black-billed magpie',
    grepl('oystercatcher', source) & grepl('Germany', location) ~ 'eurasian oystercatcher',
    grepl('canade goose', source) ~ 'canada goose',
    grepl('european herring gull', source) ~ 'herring gull',
    grepl('towny owel', source) ~ 'tawny owl',
    grepl('gadwall duck', source) ~ 'gadwall',
    grepl('lesser snow goose', source) ~ 'snow goose',
    grepl("pink footed goose", source) ~ "pink-footed goose",
    grepl("european wigeon", source) ~ "eurasian wigeon",
    source == 'sea eagle' & grepl('Norway', location)~"white-tailed eagle",
    grepl('guinea fowl', source) & grepl('Italy|Estonia|Germany', location) ~ "helmeted guineafowl (domestic type)",
    grepl('western jackdaw', source) ~ 'eurasian jackdaw',
    grepl('^gannet$', source) & grepl('Scotland', location) ~ "northern gannet",
    grepl('spotbill duck', source) ~ "indian spot-billed duck",
    grepl("whistling duck", source) ~ "whistling-duck sp.",
    grepl('pink footed goose', source) ~ 'pink-footed goose',
    grepl('falcated teal', source) ~ 'falcated duck',
    grepl('^pintail$', source) ~ 'northern pintail',
    grepl('grey plover', source)~ "black-bellied plover",
    grepl('greenwing duck', source) ~ "green-winged teal",
    grepl('franklins gull', source) ~ "franklin's gull",
    
    
    
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
    
    grepl('^seal$', source) ~ 'seal sp.', 
    grepl('^mink$', source) ~ 'mink sp.',
    grepl('^otter$', source) & grepl('Scotland', location) ~ 'eurasian otter',
    grepl('badger', source) ~ 'european badger',
    
    
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
    grepl('^pochard$', source) ~ "aythya sp.",
    grepl('coot', source) ~ 'coot sp.',
    grepl("waterfowl", source) ~ "waterfowl sp.",
    grepl("peacock|peafowl", source) ~ "indian peafowl (domestic type)",
    grepl("quail", source) & location != "America" ~ "old world quail sp.", # requires better location
    grepl("shorebird|sandpiper", source) ~ "shorebird sp.",
    grepl("black-backed gull|^gull$|seagull", source) ~ "gull sp.",
    grepl("^crow$", source) ~ "crow/raven sp.",
    grepl("^duck$|wild duck|migratory duck", source) ~ "duck sp.",
    grepl("^grebe$", source) ~ "grebe sp.",
    grepl("^owl$", source) ~ "owl sp.",
    grepl('^hawk$', source) ~ "hawk sp.",
    grepl('^eagle$', source) ~ 'eagle sp.',
    grepl('steamer duck', source) ~ "steamer-duck sp.",
    grepl('^falcon$', source) ~ "falcon sp.",
    grepl("^pelican$", source) ~ "pelican sp.",
    grepl("^pheasant$", source) ~ "pheasant sp.",
    grepl("^pigeon$", source) ~ "pigeon/dove sp.",
    grepl("^teal$", source) ~ "teal sp.",
    grepl('^ibis$', source) ~ 'ibis sp.',
    source %in% c("aquatic bird", "avian", "bird", "wild bird", "wild bird", "backyard bird") ~ "bird sp.",
    source %in% c("en", "env", "environment", "enviroment", "environmental", 'environment sample',
                  "water", "wild bird feces") ~ "environment",
    grepl("red fox|^fox$", source) ~ "red fox",
    grepl("attus norvegicus", source) ~ "brown rat",
    .default = source
  )) %>%
  
  # Get host taxonomy
  rowwise() %>%
  mutate(taxonomy = getTaxonomyForName(primary_com_name)) %>%
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

    
    # The plan is to format both frames properly, then left join to 'main' metadata, r
  