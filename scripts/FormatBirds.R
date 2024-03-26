# 
isBird <- function(x){
  bird_orders <- c(
    "passeriformes", "charadriiformes", "suliformes", "pelecaniformes","accipitriformes", 
    "falconiformes", "galliformes", "piciformes", "anseriformes", "coraciiformes", "strigiformes", 
    "columbiformes", "gruiformes","procellariiformes", "ciconiiformes", "podicipediformes", 
    "apodiformes", "caprimulgiformes", "psittaciformes", "sphenisciformes", "casuariiformes",
    "phoenicopteriformes", "cuculiformes", "coraciformes", "trogoniformes", "gaviiformes", 
    "pteroclidiformes", "anseranseriformes", "phaethontiformes","opisthocomiformes")
  
  if(x %in% bird_orders){
    x <- TRUE
  }else{
    x <- FALSE
  }
  
  return(x)
}

# input from sam's function - may be superfluous here due to domestic rows
isDomestic<- function(x){
  domestic_birds <- c(
    "ck", "ch", "chicken", "guinea_fowl", "guineafowl", "sck", "silky_chicken",
    "silkie_chicken", "silkie", "silky_fowl", "fowl", "poultry", "quail", "gs",
    "goose", "dk", "duck", "pheasant", "ph", "partridge", "gf", "turkey",
    "village_chicken", "muscovy_duck", "muscovy", "cairina_moschata", "bantam", 
    "breeder_duck", "broiler_duck", "rooster", "breeder_chicken", "mulard", 
    "mulard_duck", "korean_native_chicken", "bronze_turkey", "black_chicken", 
    "anser_cygnoides", "indian_runner_duck", "gallus_gallus",
    "meleagris_gallopavo",  "white_peacock", "peacock", "peafowl", "peahen",
    "laying_hen", "layer_hen", "chiken", 
    "broiler", "hen", "layer", "fancy_chicken", "breeding_hen", "backyard_chicken", 
    "broiler_chicken", "chichen", "layer_chicken", "rural_laying_hen", "buff_orpington_chicken", 
    "brahma_chicken", "guinea_fowl-chicken")
  
  if(x %in% domestic_birds){
    x <- TRUE
  }else{
    x <- FALSE
  }
  
  return(x)
}


FormatAnseriformes <- function(x){
  passeriformes <- c(
    "whooper swan" = "cygnus cygnus|whooper swan",
    "mute swan" = "cygnus olor|mute swan",
    "domestic goose sp. (domestic type)" = "[domestic,pomeranian,embden,sebastopol,american buff,african,rural] goose|anser anser domesticus",
    "mallard (domestic type)" = "[domestic,pekin,mule,runner,cascade,rural,mulard] duck",
    "barnacle goose" = "branta leucopsis|barnacle goose",
    "mandarin duck" = "aix galericulata|mandarin duck",
    "duck sp." = "anatidae",
    "greater/lesser white-fronted goose" = "^white[- ]fronted goose$",
    "goose sp." = "anatidae (goose sp.)",
    "greylag goose" = "anser anser|^[grey,grey {0,1}lag,gray,gray {0,1}lag] goose$|grey go*",
    "green-winged teal" = "anas crecca|greenwing duck|[common,eurasian][ -]teal|^green[- ]{0,1}winged[- ]{0,1}teal$",
    "muscovy duck" = "cairina moschata|muscovy duck", #
    "mallard" = "anas platyrhynchos|mallard duck|^mal$|plath{0,1}yrh{0,1}ynchos|^mallard$",
    "swan sp." = "^cygnus$|^swan$",
    "greater white-fronted goose" = "anser albifrons|greater[- ]white[- ]fronted[- ]goose",
    "northern pintail" = "anas acuta|^pintail$|northern pintail",
    "canada goose" = "branta canadensis|canad[ae] goose",
    "whistling-duck sp." = "dendrocygna|whistling[- ]duck",
    "falcated duck" = "mareca falcata|^falcated [teal,duck]$",
    "brant" = "branta bernicla|^brant$",
    "egyptian goose" = "alopochen aegyptiaca|egyptian goose",
    "gadwall" = "mareca strepera|^gadwall*",
    "pink-footed goose" = "anser brachyrhynchus|pink[- ]footed goose",
    "black swan" = "cygnus atratus|[black,bk][- ]swan",
    "tundra swan" = "cygnus columbianus|^[tundra,bewicks{0,1}] swan$",
    "ferruginous duck" = "aythya nyroca|ferruginous duck",
    "garganey" = "spatula querquedula|garganey",
    "taiga/tundra bean-goose" = "anser fabalis|anser serrirostris|bean[ -]goose|taiga/tundra bean-goose",
    "teal sp." = "^teal$",
    "eastern spot-billed duck" = "anas zonorhyncha|spot-billed duck|eastern spot-billed duck",
    "tufted duck" = "aythya fuligula|tufted duck",
    "eurasian/american wigeon" = "mareca [penelope,americana]|^wigeon$",
    "hawaiian goose" = "branta sandvicensis|[hawaiian,nene] goose",
    "waterfowl sp." = "anseriformes|waterfowl|wild waterbird",
    "swan goose" = "anser cygnoides|swan goose",
    "american wigeon" = "mareca americana|american wigeon",
    "northern shoveler" = "spatula clypeata|[northern ]{0,1}shoveler",
    "bufflehead" = "bucephala albeola|bufflehead",
    "wood duck" = "aix sponsa|[wood,carolina] duck",
    "lesser scaup" = "aythya affinis|lesser scaup",
    "snow goose" = "anser caerulescens|snow goose",
    "blue-winged teal" = "spatula discors|blue[ -]winged[ -]teal",
    "red-breasted goose" = "branta ruficollis|red[- ]breasted goose",
    "ross's goose" = "anser rossii|ross'{0,1}s goose",
    "merganser sp." = "mergus|merganser",
    "dabbling duck sp." = "shelduck",
    "bar-headed goose" = "anser indicus|bar[- ]headed goose",
    "cape teal" = "anas capensis|cape teal",
    "coscoroba swan" = "coscoroba coscoroba|coscoroba swan",
    "trumpeter swan" = "cygnus buccinator|trumpeter swan",
    "cackling goose" = "branta hutchinsii|cackling goose")

  for (i in 1:length(passeriformes)){
    if(any(grepl(passeriformes[[i]], x, fixed = T))){
      x <- names(passeriformes)[[i]]
    }
  }
  
  return(x)
}


FormatCharadriiformes <- function(x){
  charadriiformes <- c(
    "eurasian curlew" = "numenius arquata|eurasian curlew",
    "black-headed gull" = "chroicocephalus ridibundus|black-headed gull",
    "eurasian oystercatcher" = "haematopus ostralegus|eurasian oystercatcher",
    "lapwing sp." = "vanellus sp.|lapwing sp.",
    "red knot" = "calidris canutus|red knot",
    "herring gull" = "larus argentatus|herring gull",
    "sanderling" = "calidris alba|sanderling",
    "oystercatcher sp." = "haematopus sp.|oystercatcher sp.",
    "gull sp." = "larus sp.|gull sp.",
    "brown-headed gull" = "chroicocephalus brunnicephalus|brown-headed gull",
    "great skua" = "stercorarius skua|great skua",
    "curlew sp." = "numenius sp.|curlew sp.",
    "yellow-legged gull" = "larus michahellis|yellow-legged gull",
    "whiskered tern" = "chlidonias hybrida|whiskered tern",
    "ruddy turnstone" = "arenaria interpres|ruddy turnstone",
    "gull/tern sp." = "laridae/ternidae sp.|gull/tern sp.",
    "caspian gull" = "larus cachinnans|caspian gull",
    "common gull" = "larus canus|common gull",
    "hartlaub's gull" = "chroicocephalus hartlaubii|hartlaub's gull",
    "kelp gull" = "larus dominicanus|kelp gull",
    "long-tailed jaeger" = "stercorarius longicaudus|long-tailed jaeger",
    "little gull" = "hydrocoloeus minutus|little gull",
    "sandwich tern" = "thunnus alalunga|sandwich tern",
    "great crested tern" = "thalasseus bergii|great crested tern",
    "african oystercatcher" = "haematopus moquini|african oystercatcher",
    "arctic tern" = "sterna paradisaea|arctic tern",
    "ring-billed gull" = "larus delawarensis|ring-billed gull",
    "franklin's gull" = "larus pipixcan|franklin's gull",
    "common tern" = "sterna hirundo|common tern",
    "slaty-backed gull" = "larus schistisagus|slaty-backed gull",
    "black skimmer" = "rynchops niger|black skimmer",
    "glaucous-winged gull" = "larus glaucescens|glaucous-winged gull",
    "black-legged/red-legged kittiwake" = "rissa tridactyla|black-legged/red-legged kittiwake",
    "mediterranean gull" = "ichthyaetus melanocephalus|mediterranean gull",
    "slender-billed gull" = "larus genei|slender-billed gull",
    "glaucous gull" = "larus hyperboreus|glaucous gull",
    "western gull" = "larus occidentalis|western gull"
  )
  
  for (i in 1:length(charadriiformes)){
    if(any(grepl(charadriiformes[[i]], x, fixed = T))){
      x <- names(charadriiformes)[[i]]
    }
  }
  
  return(x)
}



FormatAccipitriformes <- function(x){
  accipitriformes <- c(
    "hawk sp." = "accipitridae sp.|hawk sp.",
    "common buzzard" = "buteo buteo|common buzzard",
    "eastern buzzard" = "buteo japonicus|eastern buzzard",
    "eurasian goshawk" = "accipiter gentilis|eurasian goshawk",
    "white-tailed eagle" = "haliaeetus albicilla|white-tailed eagle",
    "golden eagle" = "aquila chrysaetos|golden eagle",
    "western marsh harrier" = "circus aeruginosus|western marsh harrier",
    "eagle sp." = "accipitridae (eagle sp.)|eagle sp.",
    "eurasian sparrowhawk" = "accipiter nisus|eurasian sparrowhawk",
    "cooper's hawk" = "accipiter cooperii|cooper's hawk",
    "bald eagle" = "haliaeetus leucocephalus|bald eagle",
    "red-tailed hawk" = "buteo jamaicensis|red-tailed hawk",
    "african fish-eagle" = "haliaeetus vocifer|african fish-eagle",
    "old world vulture sp." = "accipitridae (old world vulture sp.)|old world vulture sp.",
    "variable hawk" = "buteo polyosoma|variable hawk",
    "bearded vulture" = "gypaetus barbatus|bearded vulture",
    "red-shouldered hawk" = "buteo lineatus|red-shouldered hawk",
    "osprey" = "pandion haliaetus|osprey",
    "harris's hawk" = "parabuteo unicinctus|harris's hawk",
    "rough-legged hawk" = "buteo lagopus|rough-legged hawk",
    "swainson's hawk" = "buteo swainsoni|swainson's hawk"
  )
  
  for (i in 1:length(accipitriformes)){
    if(any(grepl(accipitriformes[[i]], x, fixed = TRUE))){
      x <- names(accipitriformes)[[i]]
    }
  }
  
  return(x)
}

FormatPelecaniformes <- function(x){
  pelecaniformes <- c(
    "great egret" = "ardea alba|great egret",
    "grey heron" = "ardea cinerea|grey heron",
    "white egret sp." = "egretta sp.|white egret sp.",
    "great white pelican" = "pelecanus onocrotalus|great white pelican",
    "dalmatian pelican" = "pelecanus crispus|dalmatian pelican",
    "pelican sp." = "pelecanidae (pelican sp.)|pelican sp.",
    "eurasian spoonbill" = "platalea leucorodia|eurasian spoonbill",
    "ibis sp." = "threskiornithidae (ibis sp.)|ibis sp.",
    "heron sp." = "ardeidae (heron sp.)|heron sp.",
    "black-headed heron" = "ardea melanocephala|black-headed heron",
    "american white pelican" = "pelecanus erythrorhynchos|american white pelican",
    "little egret" = "egretta garzetta|little egret",
    "great blue heron" = "ardea herodias|great blue heron",
    "brown pelican" = "pelecanus occidentalis|brown pelican",
    "snowy egret" = "egretta thula|snowy egret"
  )
  
  for (i in 1:length(pelecaniformes)){
    if(any(grepl(pelecaniformes[[i]], x, fixed = TRUE))){
      x <- names(pelecaniformes)[[i]]
    }
  }
  
  return(x)
}


FormatGalliformes <- function(x){
  galliformes <-  c(
    "red junglefowl (domestic type)" = "gallus gallus|red junglefowl (domestic type)",
    "wild turkey (domestic type)" = "meleagris gallopavo|wild turkey (domestic type)",
    "helmeted guineafowl (domestic type)" = "numida meleagris|helmeted guineafowl (domestic type)",
    "pheasant sp." = "phasianidae (pheasant sp.)|pheasant sp.",
    "indian peafowl (domestic type)" = "pavo cristatus|indian peafowl (domestic type)",
    "old world quail sp." = "phasianidae (quail sp.)|old world quail sp.",
    "wild turkey" = "meleagris gallopavo|wild turkey",
    "new world quail sp." = "odontophoridae (quail sp.)|new world quail sp.",
    "lady amherst's pheasant" = "chrysolophus amherstiae|lady amherst's pheasant",
    "ring-necked pheasant" = "phasianus colchicus|ring-necked pheasant",
    "crested guineafowl sp." = "guttera pucherani|crested guineafowl sp.",
    "indian peafowl" = "pavo cristatus|indian peafowl"
  )
  
  for (i in 1:length(galliformes)){
    if(any(grepl(galliformes[[i]], x, fixed = TRUE))){
      x <- names(galliformes)[[i]]
    }
  }
  
  return(x)
}

FormatStrigiformes <- function(x){
  strigiformes <- c(
    "owl sp." = "strigiformes (owl sp.)|owl sp.",
    "short-eared owl" = "asio flammeus|short-eared owl",
    "eurasian eagle-owl" = "bubo bubo|eurasian eagle-owl",
    "tawny owl" = "strix aluco|tawny owl",
    "long-eared owl" = "asio otus|long-eared owl",
    "little owl" = "athene noctua|little owl",
    "barn owl" = "tyto alba|barn owl",
    "great horned owl" = "bubo virginianus|great horned owl",
    "ural owl" = "strix uralensis|ural owl",
    "eurasian scops-owl" = "otus scops|eurasian scops-owl",
    "snowy owl" = "bubo scandiacus|snowy owl",
    "western screech-owl" = "megascops kennicottii|western screech-owl"
  )
  
  for (i in 1:length(strigiformes)){
    if(any(grepl(strigiformes[[i]], x, fixed = TRUE))){
      x <- names(strigiformes)[[i]]
    }
  }
  
  return(x)
}

FormatPasseriformes <- function(x){
  passeriformes <-  c(
    "eurasian jackdaw" = "corvus monedula|eurasian jackdaw",
    "eurasian jay" = "garrulus glandarius|eurasian jay",
    "eurasian magpie" = "pica pica|eurasian magpie",
    "house sparrow" = "passer domesticus|house sparrow",
    "crow/raven sp." = "corvidae (crow/raven sp.)|crow/raven sp.",
    "large-billed crow" = "corvus macrorhynchos|large-billed crow",
    "song thrush" = "turdus philomelos|song thrush",
    "black-billed magpie" = "pica hudsonia|black-billed magpie",
    "american crow" = "corvus brachyrhynchos|american crow",
    "acrocephalus sp." = "acrocephalus (acrocephalus sp.)|acrocephalus sp.",
    "crow sp." = "corvus (crow sp.)|crow sp.",
    "great-tailed grackle" = "quiscalus mexicanus|great-tailed grackle"
  )
  
  for (i in 1:length(passeriformes)){
    if(any(grepl(passeriformes[[i]], x, fixed = TRUE))){
      x <- names(passeriformes)[[i]]
    }
  }
  
  return(x)
}

FormatSuliformes <- function(x){
  suliformes <- c(
    "cormorant sp." = "phalacrocoracidae (cormorant sp.)|cormorant sp.",
    "great cormorant" = "phalacrocorax carbo|great cormorant",
    "northern gannet" = "morus bassanus|northern gannet",
    "cape cormorant" = "phalacrocorax capensis|cape cormorant",
    "sulid sp." = "sulidae (sulid sp.)|sulid sp.",
    "cape gannet" = "morus capensis|cape gannet",
    "peruvian booby" = "sula variegata|peruvian booby",
    "guanay cormorant" = "leucocarbo bougainvillii|guanay cormorant",
    "neotropic cormorant" = "phalacrocorax brasilianus|neotropic cormorant",
    "double-crested cormorant" = "phalacrocorax auritus|double-crested cormorant",
    "brown booby" = "sula leucogaster|brown booby"
  )
  
  for (i in 1:length(suliformes)){
    if(any(grepl(suliformes[[i]], x, fixed = TRUE))){
      x <- names(suliformes)[[i]]
    }
  }
  
  return(x)
}

FormatGruiformes <- function(x){
  gruiformes <- c()
  
  for (i in 1:length(gruiformes)){
    if(any(grepl(gruiformes[[i]], x, fixed = TRUE))){
      x <- names(gruiformes)[[i]]
    }
  }
  
  return(x)
}

FormatPodicipediformes <- function(x){
  podicipediformes <- c()
  
  for (i in 1:length(podicipediformes)){
    if(any(grepl(podicipediformes[[i]], x, fixed = TRUE))){
      x <- names(podicipediformes)[[i]]
    }
  }
  
  return(x)
}

FormatFalconiformes <- function(x){
  falconiformes <- c()
  
  for (i in 1:length(falconiformes)){
    if(any(grepl(falconiformes[[i]], x, fixed = TRUE))){
      x <- names(falconiformes)[[i]]
    }
  }
  
  return(x)
}

FormatCiconiiformes <- function(x){
  ciconiiformes <- c()
  
  for (i in 1:length(ciconiiformes)){
    if(any(grepl(ciconiiformes[[i]], x, fixed = TRUE))){
      x <- names(ciconiiformes)[[i]]
    }
  }
  
  return(x)
}

FormatRheiformes <- function(x){
  rheiformes <- c()
  
  for (i in 1:length(rheiformes)){
    if(any(grepl(rheiformes[[i]], x, fixed = TRUE))){
      x <- names(rheiformes)[[i]]
    }
  }
  
  return(x)
}

FormatColumbiformes <- function(x){
  columbiformes <- c()
  
  for (i in 1:length(columbiformes)){
    if(any(grepl(columbiformes[[i]], x, fixed = TRUE))){
      x <- names(columbiformes)[[i]]
    }
  }
  
  return(x)
}

FormatPsittaciformes <- function(x){
  psittaciformes <- c()
  
  for (i in 1:length(psittaciformes)){
    if(any(grepl(psittaciformes[[i]], x, fixed = TRUE))){
      x <- names(psittaciformes)[[i]]
    }
  }
  
  return(x)
}

FormatCathartiformes <- function(x){
  cathartiformes <- c()
  
  for (i in 1:length(cathartiformes)){
    if(any(grepl(cathartiformes[[i]], x, fixed = TRUE))){
      x <- names(cathartiformes)[[i]]
    }
  }
  
  return(x)
}


FormatProcellariiformes <- function(x){
  procellariiformes <- c()
  
  for (i in 1:length(procellariiformes)){
    if(any(grepl(procellariiformes[[i]], x, fixed = TRUE))){
      x <- names(procellariiformes)[[i]]
    }
  }
  
  return(x)
}

FormatSphenisciformes <- function(x){
  sphenisciformes <- c()
  
  for (i in 1:length(sphenisciformes)){
    if(any(grepl(sphenisciformes[[i]], x, fixed = TRUE))){
      x <- names(sphenisciformes)[[i]]
    }
  }
  
  return(x)
}

FormatCasuariiformes <- function(x){
  casuariiformes <- c()
  
  for (i in 1:length(casuariiformes)){
    if(any(grepl(casuariiformes[[i]], x, fixed = TRUE))){
      x <- names(casuariiformes)[[i]]
    }
  }
  
  return(x)
}

FormatStruthioniformes <- function(x){
  struthioniformes <- c()
  
  for (i in 1:length(struthioniformes)){
    if(any(grepl(struthioniformes[[i]], x, fixed = TRUE))){
      x <- names(struthioniformes)[[i]]
    }
  }
  
  return(x)
}

FormatPhoenicopteriformes <- function(x){
  phoenicopteriformes <- c()
  
  for (i in 1:length(phoenicopteriformes)){
    if(any(grepl(phoenicopteriformes[[i]], x, fixed = TRUE))){
      x <- names(phoenicopteriformes)[[i]]
    }
  }
  
  return(x)
}



#FormatCoraciformes
#FormatBucerotiformes
#FormatTrogoniformes
#FormatLeptosomiformes
#FormatColilformes
#FormatOpisthocomiformes
#FormatGaviiformes
#FormatPterocliformes
#FormatCucliformes
#FormatOtidiformes
#FormatMusophagiformes
#FormatApodiformes
#FormatCaprimulgiformes
#FormatApterygiformes
#FormatTinamiformes
