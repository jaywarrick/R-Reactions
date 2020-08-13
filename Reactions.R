Reaction <- function(name)
{
    reaction <- new('Reaction', name=name)
    setParams(reaction, params=paste('k_',name,sep=''))
    return(reaction)
}


setClass('Reaction', 
         slots = c(name='character', reactants='data.frame', products='data.frame', params='character', baserate.eq='character', rates='data.frame'),
         # prototype = c(
         #     name='',
         #     reactants=data.frame(name=character(0), stoich=numeric(0), order=numeric(0), stringsAsFactors=FALSE),
         #     products=data.frame(name=character(0), stoich=numeric(0), stringsAsFactors=FALSE),
         #     params=character(0),
         #     baserate.eq='',
         #     rates=data.frame(name=character(0), stringsAsFactors=FALSE)
         # )
)

setGeneric('setParams',function(object='Reaction', params='character'){standardGeneric ('setParams')})
setMethod(f = 'setParams', signature = c(object='Reaction', params='character'), 
          definition = function(object, params)
          {
              passedObject<-deparse(substitute(object))
              object@params <- params
              assign(passedObject,object,envir=parent.frame())
          })

setGeneric('setReactant',function(object='Reaction', name='character', stoich='numeric', order='numeric'){standardGeneric ('setReactant')})
setMethod(f = 'setReactant', signature = c(object='Reaction', name='character', stoich='numeric', order='numeric'), 
          definition = function(object, name, stoich=1, order=1)
          {
              passedObject<-deparse(substitute(object))
              if(length(which(object@reactants$name==name) > 0))
              {
                  object@reactants$stoich <- stoich
                  object@reactants$order <- order
              }
              else
              {
                  object@reactants <- rbind(object@reactants, data.frame(name=name, stoich=stoich, order=order, stringsAsFactors=FALSE))
              }
              assign(passedObject,object,envir=parent.frame())
          })

setGeneric('setProduct',function(object='Reaction', name='character', stoich='numeric'){standardGeneric ('setProduct')})
setMethod(f = 'setProduct', signature(object='Reaction', name='character', stoich='numeric'), 
          definition = function(object, name, stoich=1)
          {
              passedObject<-deparse(substitute(object))
              if(length(which(object@products$name==name) > 0))
              {
                  object@products$stoich <- stoich
              }
              else
              {
                  object@products <- rbind(object@products, data.frame(name=name, stoich=stoich, stringsAsFactors=FALSE))
              }
              assign(passedObject,object,envir=parent.frame())
          })

setGeneric('setBaseRateEq',function(object='Reaction', eq='character'){standardGeneric ('setBaseRateEq')})
setMethod(f = 'setBaseRateEq', signature(object='Reaction', eq='character'),
          definition = function(object, eq='0')
          {
              passedObject<-deparse(substitute(object))
              object@baserate.eq <- eq
              assign(passedObject,object,envir=parent.frame())
          })

setGeneric('makeRates',function(object='Reaction'){standardGeneric ('makeRates')})
setMethod(f = 'makeRates', signature(object='Reaction'), 
          definition = function(object)
          {
              passedObject<-deparse(substitute(object))
              if(object@baserate.eq=='')
              {
                  baserate <- object@params[1]
                  for(r in object@reactants$name)
                  {
                      reactant <- subset(object@reactants, name==r)
                      baserate <- paste(baserate, '*', reactant$name, '^(', reactant$order, ')', sep='')
                  }
                  object@baserate.eq <- baserate
              }
              rates <- data.frame(name=character(0), eq=character(0), stringsAsFactors=FALSE)
              baserate <- object@baserate.eq
              cancelations <- character(0);
              for(r in object@reactants$name)
              {
                  reactant <- subset(object@reactants, name==r)
                  matchingIndex <- which(object@products$name==r)                  
                  hasMatch <- !(length(matchingIndex)==0)
                  stoichMatches <- FALSE
                  if(hasMatch)
                  {
                      stoichMatches <- (object@products$stoich[matchingIndex]==reactant$stoich)
                  }
                  if(!hasMatch | !stoichMatches)
                  { # then add the reactant consumption because it isn't added back equally as a product
                      eq <- paste('-',reactant$stoich, '*', baserate, sep='')
                      name <- paste('d', reactant$name, sep='')
                      rates <- rbind(rates, data.frame(name=name, eq=eq, stringsAsFactors=FALSE))
                  }
                  else
                  {
                      # the reactant shows up as a reactant and product in equal stoichiometry and should be omitted
                      cancelations <- c(cancelations, reactant$name)
                  }
              }
              for(p in object@products$name)
              {
                  product <- subset(object@products, name==p)
                  cancelation <- length(which(cancelations == p)) > 0
                  if(!cancelation)
                  {
                      eq <- paste(product$stoich, '*', baserate, sep='')
                      name <- paste('d', product$name, sep='')
                      rates <- rbind(rates, data.frame(name=name, eq=eq, stringsAsFactors=FALSE))
                  }
                  else
                  {
                      eq <- '0'
                      name <- paste('d', product$name, sep='')
                      rates <- rbind(rates, data.frame(name=name, eq=eq, stringsAsFactors=FALSE))
                  }
              }
              object@rates <- rates
              assign(passedObject,object,envir=parent.frame())
          })


setClass('ODEModel', 
         slots = c(name='character', eqs='data.frame', vars='character', preCalcCode='character'),
         # prototype = c(
         #     name='',
         #     eqs=data.frame(name=character(0), eq=character(0), stringsAsFactors=FALSE),
         #     vars=character(0),
         #     preCalcCode=''
         # )
)

ODEModel <- function(name)
{
    ret <- new('ODEModel')
    return(ret)
}

setGeneric('addReaction',function(object='ODEModel', reaction='Reaction'){standardGeneric ('addReaction')})
setMethod(f = 'addReaction', signature(object='ODEModel', reaction='Reaction'),
          definition = function(object, reaction)
          {
              passedObject<-deparse(substitute(object))
              eqs <- object@eqs
              
              makeRates(reaction)
              theRates <- reaction@rates
              alreadyHave <- FALSE
              for(i in 1:nrow(theRates))
              {
                  theName <- theRates$name[i]
                  if(nrow(eqs) == 0)
                  {
                      alreadyHave <- FALSE
                  }
                  else
                  {
                      alreadyHave <- (length(which(eqs$name==theName)) > 0)
                  }
                  
                  if(alreadyHave)
                  {
                      eqs[eqs$name==theName,]$eq <- paste(eqs$eq[eqs$name==theName], ' + ', theRates$eq[i], sep='')
                  }
                  else
                  {
                      eqs <- rbind(eqs, theRates[i,])
                  }
              }
              
              object@eqs <- eqs
              assign(passedObject,object,envir=parent.frame())
          })

setGeneric('addDE',function(object='ODEModel', dX='character', eq='character'){standardGeneric ('addDE')})
setMethod(f = 'addDE', signature(object='ODEModel', dX='character', eq='character'),
          definition = function(object, dX, eq)
          {
              passedObject<-deparse(substitute(object))
              eqs <- object@eqs
              
              alreadyHave <- FALSE
              theName <- dX
              if(nrow(eqs) == 0)
              {
                  alreadyHave <- FALSE
              }
              else
              {
                  alreadyHave <- (length(which(eqs$name==theName)) > 0)
              }
              
              if(alreadyHave)
              {
                  eqs[eqs$name==theName,]$eq <- paste(eqs$eq[eqs$name==theName], ' + ', eq, sep='')
              }
              else
              {
                  eqs <- rbind(eqs, data.frame(name=theName, eq=eq, stringsAsFactors=FALSE))
              }
              
              object@eqs <- eqs
              assign(passedObject,object,envir=parent.frame())
          })

setGeneric('setPreCalcCode',function(object='ODEModel', code='character'){standardGeneric ('setPreCalcCode')})
setMethod(f = 'setPreCalcCode', signature(object='ODEModel', code='character'),
          definition = function(object, code)
          {
               passedObject<-deparse(substitute(object))
               object@preCalcCode <- code
               assign(passedObject,object,envir=parent.frame())
          })

setGeneric('makeODEFunc',function(object='ODEModel'){standardGeneric ('makeODEFunc')})
setMethod(f = 'makeODEFunc', signature(object='ODEModel'),
          definition = function(object)
          {
              passedObject<-deparse(substitute(object))
               
              # Assemble the header
              header <- 'function(t, X, params) {\n with(as.list(c(X, params)), {\n'
              
              # Assemble pre-calc equations
              prebody <- object@preCalcCode
              
              # Assemble the DE rate equations
              body <- ''
              object@vars <- character(0)
              eqs <- object@eqs
              for(item in eqs$name)
              {
                  body <- paste(body, item, ' = ', subset(eqs, name==item)$eq, '\n', sep='')
                  object@vars <- c(object@vars, substr(item, 2, nchar(item)))
              }
              dVars <- sapply(object@vars,FUN=function(x){paste('d',x,sep='')})
              returnStatement <- paste('return(list(c(', paste(dVars, collapse=','), ')))')
              close <- '})}'
              
              # Assemble the text of the function and compile it
              print(eqs)
              theFunction <- paste(header, '\n', prebody, '\n', body, '\n', returnStatement, '\n', close, sep='')
              assign(passedObject,object,envir=parent.frame())
              eval(parse(text=theFunction))
          })

setGeneric('getOrderedICVector',function(object='ODEModel', ic='list'){standardGeneric ('getOrderedICVector')})
setMethod(f = 'getOrderedICVector', signature(object='ODEModel', ic='list'),
          definition = function(object, ic)
          {
               passedObject<-deparse(substitute(object))
               orderedICs <- list()
               vars <- object@vars
               for(varName in vars)
               {
                    orderedICs[varName] <- ic[varName]
               }
               print(unlist(orderedICs))
               return(unlist(orderedICs))
          })







