# Load required libraries
import synapseclient
import synapseutils

# login to Synapse
syn = synapseclient.login(
  email='', # your synapse email id
  password='' # your password
) 
files = synapseutils.syncFromSynapse(syn, entity = 'syn20564743', path = './data/')

