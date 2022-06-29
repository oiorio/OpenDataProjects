import FWCore.ParameterSet.Config as cms

slistttdl= cms.untracked.vstring()

slistttdl.extend(
["root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/00610F90-B8EC-E411-8BB0-0025905A48F0.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/042E49D5-A4EC-E411-BB3D-0025905A48E4.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0A2B1CCA-4DEC-E411-9A0B-0025905B8572.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0AB4F55B-A4EC-E411-9D9C-0025905B85AA.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0C0844E9-0EEC-E411-A191-002590596468.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/0C462BC1-4DEC-E411-9475-003048FFCC1E.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/108C7D61-F7EC-E411-9CDB-0025905A6064.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/14C7580F-ADEC-E411-9607-0025905B8572.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/18951A1B-9AEC-E411-A7C6-00261894393E.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/1A86A528-0AEC-E411-BB58-0025905A6090.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/1AED8B9B-ADEC-E411-8D4B-003048FFD7D4.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/1E54947F-ABEC-E411-A831-0025905A48F2.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/28911B56-4DEC-E411-9B85-003048FFCC18.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/2A5932A7-70EC-E411-AF35-00259059642A.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/2C61CD5B-4EEC-E411-81C2-002618943981.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/2CC6D6E1-AEEC-E411-A967-0025905B85AE.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/30AAF936-AEEC-E411-B16E-0025905A48EC.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/3228430F-ADEC-E411-94DC-0025905A612E.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/32CD9586-ACEC-E411-B616-0025905B85EE.root",
"root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2011/Summer11LegDR/TTJets_MSDecays_dileptonic_central_7TeV-madgraph-tauola/AODSIM/PU_S13_START53_LV6-v1/00000/364A12ED-16EC-E411-8848-0025905B858E.root",])
