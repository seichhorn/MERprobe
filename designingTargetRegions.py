import sys
import pandas as pd
from collections import defaultdict
import pdb
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
import networkx as nx


class transcript:
  """A class to contain information about individual transcripts"""
  def __init__(self, identifier,sequence):
    self.identifier = identifier
    self.sequence = sequence
    self.expression = 0.0
    self.geneName = identifier

  def updateExpression(self,expression):
    self.expression = expression

  def updateGeneName(self,geneName):
    self.geneName = geneName

class transcriptome:
  """A class to contain a collection of transcripts"""
  def __init__(self):
    self.transcripts = dict()
    self.numberOfTranscripts = len(self.transcripts)

  def addTranscript(self,transcript):
    if transcript in self.transcripts:
      print('Attempted to add a transcript that is already present')
    else:
      self.transcripts[transcript.identifier] = transcript
      self.numberOfTranscripts = len(self.transcripts)

  def listContainedTranscripts(self):
    self.geneList = list(self.transcripts.keys())
  
  def checkTranscript(self,identifier):
    if not len(self.geneList) == len(self.transcripts):
      listContainedTranscripts()
    if identifier in self.geneList:
      print(str(identifier)+' is present in the transcriptome')
    else:
      print(str(identifier)+' is absent from the transcriptome')

  def addExpressionValues(self,expressionFile):
    """Right now this expects a file that is just Transcript_ID \t expression \n, no header"""
    with open(expressionFile,'r') as expressionFileOpen:
      for line in expressionFileOpen:
        linesplit = line.rstrip('\n').split('\t')
        if len(linesplit) != 2:
          print('File is not in the expected format, please provide a 2-column file, first column is identifier, second is expression')
          break
        else:
          identifier = linesplit[0]
          expression = linesplit[1]
          if identifier in self.transcripts:
            self.transcripts[identifier].updateExpression(expression)

  def addGeneNames(self,xrefFile):
    """Right now this expects a file that is just Transcript_ID \t Gene_name \n with no header"""
    with open(xrefFile,'r') as xrefFileOpen:
      for line in xrefFileOpen:
        linesplit = line.rstrip('\n').split('\t')
        if len(linesplit) != 2:
          print('File is not in the expected format, please provide a 2-column file, first column is identifier, second is gene name')
          break
        else:
          identifier = linesplit[0]
          geneName = linesplit[1]
          if identifier in self.transcripts:
            self.transcripts[identifier].updateGeneName(geneName)

  def getTranscripts(self,listOfTranscripts):
    txCollection = []
    for tx in listOfTranscripts:
      if tx in self.transcripts:
        txCollection.append(self.transcripts[tx])
      else:
        print('Missing {} from transcriptome'.format(tx))
    return txCollection

class gene:
  def __init__(self,transcript):
    if transcript.geneName != transcript.identifier:
      self.geneName = transcript.geneName
      self.identifiers = [transcript.identifier]
      self.expression = [transcript.expression]
      self.sequences = [transcript.sequence]
      self.isoformNumber = len(self.identifiers)
    else:
      print('{} doesn\'t have a gene name'.format(transcript.identifier))
  
  def addIsoform(self,transcript):
    if self.geneName == transcript.geneName:
      self.identifiers.append(transcript.identifier)
      self.expression.append(transcript.expression)
      self.sequences.append(transcript.sequence)
      self.isoformNumber = len(self.identifiers)

class geneCollection:
  '''like a transcriptome but instead with gene-centric referencing and supports multiple isoforms'''
  def __init__(self):
    self.genes = dict()

  def assembleFromTranscriptome(self,transcriptome):
    for transcript in transcriptome.transcripts.values():
      if transcript.geneName in self.genes:
        self.genes[transcript.geneName].addIsoform(transcript)
      else:
        self.genes[transcript.geneName] = gene(transcript)

  def addGenes(self,genes):
    #genes is expected to be a list of gene class objects
    for g in genes:
      if g.geneName in self.genes:
        print('{} is already present in gene collection'.format(gene.geneName))
      else:
        self.genes[g.geneName] = g


class probe:
  """Contains the general features of a probe"""
  def __init__(self,seq,position,specificity = None):
    self.seq = seq.upper()
    self.position = position
    self.length = len(seq)
    self.gcFrac = (seq.upper().count('G') + seq.upper().count('C')) / len(seq)
    if 'X' in self.seq:
      self.tm = 0
    else:
      self.tm = mt.Tm_NN(Seq(seq),Na = 300, dnac1 = 5, dnac2 = 1, saltcorr = 7)
    if specificity == None:
      self.specificity = None
    else:
      self.specificity = specificity


class forbiddenSeq:
  """A class to contain a list of kmers that should not appear in any selected probe"""
  def __init__(self,k,chosenTranscriptome):
    self.k = k
    self.transcriptome = chosenTranscriptome
    self.forbidden = createAllPossibleTranscriptomeKmer(self.k,self.transcriptome)

def constructTranscriptomeFastA(fastAFile):
  """Creates a transcriptome based on a fastA file"""
  localTranscriptome = transcriptome()
  lastline = False
  a = 0
  b = 2
  fastAOpen = open(fastAFile,'r')
  while lastline == False:
    package = []
    for i in range(a,b):
      line = fastAOpen.readline()
      if line == '':
        lastline = True
      else:
        package.append(line)
    a += 2
    b += 2
    if len(package) == 2:
      currentGene = transcript(package[0].split('>')[1].rstrip('\n'),package[1].rstrip('\n'))
      localTranscriptome.addTranscript(currentGene)
  return localTranscriptome

def createAllPossibleKmer(k,seq):
  """Divides a sequence into all possible kmers"""
  seqSize = len(seq)
  currentPos = 0
  kmers = []
  while currentPos + k + 1 <= seqSize:
    kmers.append(seq[currentPos:currentPos+k].upper())
    currentPos+= 1
  return list(kmers)

def createAllPossibleTranscriptomeKmer(k,transcriptome):
  """Runs create all possible kmers on each element of a transcriptome"""
  kmerSet = set([])
  counter = 0
  print('{} sequences to process'.format(len(transcriptome.transcripts.values())))
  for trans in transcriptome.transcripts.values():
    counter += 1
    listOfKmers = createAllPossibleKmer(k,trans.sequence)
    kmerSet.update(set(listOfKmers))
    if counter%100 == 0:
      print('Finished {} sequences'.format(counter))
  return list(kmerSet)

def createAllPossibleProbes(k,seq):
  """Divides a sequence into all possible kmers"""
  seqSize = len(seq)
  currentPos = 0
  kmers = []
  while currentPos + k + 1 <= seqSize:
    kmers.append(probe(seq[currentPos:currentPos+k],currentPos))
    currentPos+= 1
  return list(kmers)


def makeAttributeDict(probe):
  dictionary = dict()
  dictionary['position'] = probe.position
  dictionary['length'] = probe.length
  dictionary['gcFrac'] = probe.gcFrac
  dictionary['tm'] = probe.tm
  dictionary['specificity'] = 0
  return dictionary
  

def constructGraph(transcriptList,transcriptome,k):
  """Creates a graph-based representation of all genes with edges to all possible probes"""
  g = nx.Graph()
  transcriptCollection = transcriptome.getTranscripts(transcriptList)
  for tx in transcriptCollection:
      identifier = tx.identifier
      sequence = tx.sequence
      expression = tx.expression
      kmers = createAllPossibleProbes(k,sequence)
      if len(kmers) > 0:
        attributeDict = [makeAttributeDict(x) for x in kmers]
        attributeDict = [{x['expression'] : expression} for x in attributeDict]
      edges = list(zip([identifier]*len(kmers), kmers, attributeDict))
      g.add_edges_from(edges)
  return g

def addSpecificityForSelectGenes(geneList,geneCollection,graphProbe,graphSpecificity):
  """Creates a graph-based representation of all kmers with edges to all genes"""
  geneProbeDict = defaultdict(list)
  for gene in geneList:
    for currentPosition in range(len(geneCollection.genes[gene].identifiers)):
      identifier = geneCollection.genes[gene].identifiers[currentPosition]
      expression = geneCollection.genes[gene].expression[currentPosition]
      if identifier in graphProbe:
        possibleProbes = [x[1] for x in graphProbe.edges(gene)]
        if len(possibleProbes) > 0:
          for p in possibleProbes:
            breakdown = createAllPossibleKmer(17,p)
            matchingGenes = graphSpecificity.edges(breakdown, data = 'expression')
            geneIdentifiers = geneCollection.genes[gene].identifiers
            offTarget = sum([x[2] for x in matchingGenes if x[1] not in geneIdentifiers])
            onTarget = expression
            geneProbeDict[identifier].append(probe(p,graphProbe.get_edge_data(p,identifier)['position'],onTarget/(offTarget+onTarget)))
  return geneProbeDict

def getProbesBasedOnCharacteristics(probeDict, gc = [0.4,0.6], tm = [65,75], specificity = [0.75, 1.0]):
  cutDict = defaultdict(list)
  for k,v in probeDict.items():
    for element in v:
      if element.gcFrac>=gc[0] and element.gcFrac<=gc[1] and element.tm>=tm[0] and element.tm<=tm[1]:
        cutDict[k].append(element)
  return cutDict

def selectProbesGenecentric(geneList,geneCollection,probedict):
  for gene in geneList:
    identifiers = geneCollection.genes[gene].identifiers
    expression = geneCollection.genes[gene].expression
    totalExpression = sum(expression)
    sortedIdentfier = list(sorted(zip(expression,identifiers),reverse = True))
    runningTotal = 0
    i = 0
    keepers = []
    for element in sortedIdentifiers:
      identifier = element[1]
      if i == 0:
        if identifier in cutDict:
          keepers = [(x.seq, x.position, x.specificity) for x in probedict[identifier]]
        runningTotal += element[0]
        if runningTotal/totalExpression < 0.8:
          break
      else:
        if identifier in cutDict:
          newkeepers = [(x.seq, x.position, x.specificity) for x in probedict[identifier]]
          updatedkeepers = []
          for entry in keepers:
            if entry(0) in [x[0] for x in newkeepers]:
              updatedkeepers.append(entry)
          keepers = updatedkeepers
          runningTotal += element[0]
          if runningTotal/totalExpression < 0.8:
            break

