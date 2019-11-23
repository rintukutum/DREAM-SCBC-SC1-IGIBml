# Load required libraries
import synapseclient
import synapseutils

# login to Synapse
syn = synapseclient.login(
  email='rintukutum@gmail.com',
  password='gp41@Synapse'
) 
files = synapseutils.syncFromSynapse(syn, entity = 'syn20564743', path = './data/')

