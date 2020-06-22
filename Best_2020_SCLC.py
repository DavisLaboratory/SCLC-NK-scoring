import CCLETools
import GeneSetScoring
import ENSEMBLTools
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy
from scipy.stats import gaussian_kde
import scipy.cluster.hierarchy as SciPyClus
import seaborn as sns


class PathDir:
    # A class to create folders and control folder paths

    # If desired, set this string to a path for the preproc data/output figures etc
    # strCurrDir = 'undefined'
    strCurrDir = 'C:\\git\\SCLC-NK-scoring'
    # if left undefined, this script will default to using the current working directory
    if strCurrDir == 'undefined':
        strCurrDir = os.getcwd()

    # if you wish to create specific folders for the input data, preprocesed data or figures then
    #  edit the following lines, otherwise these folders will automatically be created
    strDataLoc = os.path.join(strCurrDir, 'data')
    strTempDataLoc = os.path.join(strCurrDir, 'preproc')
    strOutputLoc = os.path.join(strCurrDir, 'figures')

    for path in [strDataLoc, strTempDataLoc, strOutputLoc]:
        if not os.path.exists(path):
            os.mkdir(path)


class Load:

    def george_2015_metadata(
            flagResult=False,
            strFileName='NIHMS782739-supplement-Tables.xlsx',
            flagPerformExtraction=False,
            strTempFileName='george_2015_metadata.pickle'):

        # This function loads the Supplementary Table (Excel spreadsheet) from George et al
        #  and exports the patient metadata sheet. As this spreadsheet can take some time to
        #  load, the script attempts to find/save a pre-processed file for improved runtime.
        #  Users can delete this temporary file or set flagPerformExtraction=True.

        # check if the pre-processed data file exists
        if not os.path.exists(os.path.join(PathDir.strTempDataLoc, strTempFileName)):
            flagPerformExtraction = True

        if flagPerformExtraction:
            print('Loading the patient metadata from George et al. (2015)')
            # Load the Excel spreadsheet and extract the corresponding sheet name
            dfPatientData = pd.read_excel(
                os.path.join(PathDir.strDataLoc, strFileName),
                sheet_name='Suppl_Table1', skiprows=3,
                header=0, index_col=0)
            # Save the temporary file as a pickle
            dfPatientData.to_pickle(
                os.path.join(PathDir.strTempDataLoc, strTempFileName))

        else:
            # Load the pre-processed pickle file
            dfPatientData = pd.read_pickle(
                os.path.join(PathDir.strTempDataLoc, strTempFileName))

        return dfPatientData


    def george_2015_transcript(
            flagResult=False,
            strFileName='NIHMS782739-supplement-Tables.xlsx',
            flagPerformExtraction=False,
            strTempFileName='george_2015_transcript.pickle'):

        # This function loads the Supplementary Table (Excel spreadsheet) from George et al
        #  and exports the transcript abundance data sheet. As this spreadsheet can take some time to
        #  load, the script attempts to find/save a pre-processed file for improved runtime.
        #  Users can delete this temporary file or set flagPerformExtraction=True.

        # check if the pre-processed data file exists
        if not os.path.exists(os.path.join(PathDir.strTempDataLoc, strTempFileName)):
            flagPerformExtraction = True

        if flagPerformExtraction:
            print('Loading the patient transcriptomic data from George et al. (2015)')
            # dfAbundData = pd.read_table(os.path.join(strDataLoc, strTranscDataFileName),
            #                             header=0, index_col=1)
            dfAbundData = pd.read_excel(
                os.path.join(PathDir.strDataLoc, strFileName),
                sheet_name='Suppl_Table10', skiprows=2,
                header=0, index_col=1)
            # Save the temporary file as a pickle
            dfAbundData.to_pickle(
                os.path.join(PathDir.strTempDataLoc, strTempFileName))

        else:
            print('Loading pre-processed patient transcriptomic from George et al. (2015)')
            # Load the pre-processed pickle file
            dfAbundData = pd.read_pickle(
                os.path.join(PathDir.strTempDataLoc, strTempFileName))

        return dfAbundData


class PreProcess:

    dictTFConstThresh = {'YAP1': 2.5,
                          'POU2F3': 2.5}

    dictTFRelBaseThresh = {'ASCL1':4.0,
                           'NEUROD1':5.0}

    dictTFRelDiffThresh = {'ASCL1':1.0,
                           'NEUROD1':1.0}


    def george_2015_small_cell(flagResult=False,
                               flagPerformExtraction=False,
                               strTempFileName='merged_data.pickle'):

        if not os.path.exists(os.path.join(PathDir.strTempDataLoc, strTempFileName)):
            flagPerformExtraction = True

        if flagPerformExtraction:
            dfPatientData = Load.george_2015_metadata()
            dfAbundData = Load.george_2015_transcript()

            print('Merging and cleaning the data from George et al. (2015) - this may take a while')
            listMetaDataPatients = \
                dfPatientData[dfPatientData['Status (at time of last follow-up)'].notnull()].index.tolist()
            listAbundPatients = dfAbundData.columns.tolist()[1:]
            print('\t .. ' + '{}'.format(len(listAbundPatients)) + ' of ' +
                  '{}'.format(len(listMetaDataPatients)) + ' patients have mRNA abundance data')

            # several patients have RNA-seq data but no available patient metadata, these include:
            #   'S01248', 'S01297', 'S01453', 'S01556', 'S01563', 'S02163', 'S02288', 'S02298'

            listPatients = listAbundPatients
            numPatients = len(listPatients)
            print('\t .. ' + '{}'.format(numPatients) + ' have valid survival data and mRNA abundance data')

            listUniqueGenes = dfAbundData['gene'].unique().tolist()

            dfData = pd.DataFrame(index=listPatients,
                                  columns=listUniqueGenes,
                                  data=np.zeros((numPatients,len(listUniqueGenes))))
            setDataGenes = set(listUniqueGenes)

            arrayGeneIsDupEntry = dfAbundData['gene'].duplicated().values.astype(np.bool)

            arrayNonUniqueRowIndices = np.where(arrayGeneIsDupEntry)[0]
            arrayUniqueRowIndices = np.where(~arrayGeneIsDupEntry)[0]

            listDupGenes = dfAbundData['gene'].iloc[arrayNonUniqueRowIndices].unique()
            listNonDupGenes = dfAbundData['gene'].iloc[arrayUniqueRowIndices].unique()

            for strGene in listNonDupGenes:
                dfData[strGene] = dfAbundData[listPatients][dfAbundData['gene'] == strGene].values.astype(np.float).transpose()

            print('\t ..  stepping through genes with duplicate entries and taking the median value')
            for strGene in listDupGenes:
                arrayAbundDataForGene = dfAbundData[listPatients][dfAbundData['gene'] == strGene]
                dfData[strGene] = np.median(arrayAbundDataForGene, axis=0)

            print('Appending metadata..')
            dfData['_smoker_status'] = dfPatientData['smoking_status'].reindex(listPatients)
            dfData['_status'] = dfPatientData['Status (at time of last follow-up)'].reindex(listPatients)
            dfData['_OS'] = dfPatientData['overall_survival (months)'].reindex(listPatients)
            dfData['_PFS'] = dfPatientData['progression-free_survival (months)'].reindex(listPatients)
            dfData['_type'] = dfPatientData['primary tumor/metastasis'].reindex(listPatients)
            dfData['_stage'] = dfPatientData['stage_UICC'].reindex(listPatients)
            dfData['_age'] = dfPatientData['age'].reindex(listPatients)
            dfData['_sex'] = dfPatientData['sex'].reindex(listPatients)

            # Load the NK cell transcriptomic signatures
            dictLocalNKSig = GeneSetScoring.CuratedList.nk_cells_in_tumour()
            listNKUpGenesInData = list(set(dictLocalNKSig['UpGenes']).intersection(setDataGenes))

            dictCombLM22TCellSig = GeneSetScoring.CuratedList.t_cells_in_tumour()
            listTCellUpGenesInData = list(set(dictCombLM22TCellSig['T cells']['UpGenes']).intersection(setDataGenes))

            print('Scoring patient samples for specified signatures..')
            arrayNKScore = np.zeros(numPatients, dtype=np.float)
            arrayTCellScore = np.zeros(numPatients, dtype=np.float)


            for iPatient in range(numPatients):
                print('{}'.format(iPatient+1) + ' of ' + '{}'.format(numPatients+1))
                arrayNKScore[iPatient] = \
                    GeneSetScoring.FromInput.single_sample_rank_score(
                        listAllGenes=listUniqueGenes,
                        arrayTranscriptAbundance=dfData[listUniqueGenes].iloc[iPatient].values.astype(np.float),
                        listUpGenesToScore=listNKUpGenesInData,
                        flagApplyNorm=True)
                arrayTCellScore[iPatient] = \
                    GeneSetScoring.FromInput.single_sample_rank_score(
                        listAllGenes=listUniqueGenes,
                        arrayTranscriptAbundance=dfData[listUniqueGenes].iloc[iPatient].values.astype(np.float),
                        listUpGenesToScore=listTCellUpGenesInData,
                        flagApplyNorm=True)

            dfData['NK score'] = arrayNKScore
            dfData['T Cell score'] = arrayTCellScore

            dfData.to_pickle(os.path.join(PathDir.strTempDataLoc, strTempFileName))

        else:

            print('Loading the pre-processed George et al. (2015) data for small cell lung cancer')
            dfData = pd.read_pickle(os.path.join(PathDir.strTempDataLoc, strTempFileName))

        return dfData


    def mollaoglu_rnaseq(flagResult=False,
                         strDataFileName='GSE89660_RPM_RPR2_fpkm.txt'):

        # Load in the original processed data frame (FPKM data from Mollaoglu et al)
        dfIn = pd.read_table(os.path.join(PathDir.strDataLoc, strDataFileName),
                             sep='\t', header=0, index_col=0)

        arrayNumObsAboveAbundThresh = np.sum(np.nan_to_num(dfIn.values.astype(np.float)) >= 0.01, axis=1)

        arrayKeepGene = arrayNumObsAboveAbundThresh >= 3

        dfData = dfIn.iloc[np.where(arrayKeepGene)[0]].copy()

        return dfData

    def mollaoglu_microarray(
            flagResult=False,
            strDataFileName='expression.txt',
            strAnnotFileName='annotation.txt'):

        dfAnnot = pd.read_table(os.path.join(PathDir.strDataLoc, strAnnotFileName),
                                sep='\t', quoting=2, header=0, index_col=None)
        listUniqueMsENSEMBLGenes = sorted(dfAnnot['ENSEMBL'].unique().tolist())
        numUniqueGenes = len(listUniqueMsENSEMBLGenes)

        # Load in the original processed data frame (FPKM data from Mollaoglu et al)
        dfIn = pd.read_table(os.path.join(PathDir.strDataLoc, strDataFileName),
                             sep='\t', header=0, index_col=0)
        listDataColumns = dfIn.columns.tolist()
        numMice = len(listDataColumns)

        listDataColumnsClean = []
        for strCol in listDataColumns:
            if 'Kwon_Sik_Park_KPPR_mSCLC' in strCol:
                strType='TKO'
            elif 'Kwon_Sik_Park_KPPR_mNomal' in strCol:
                strType='Normal'
            else:
                strType='DKO'
            strColClean = strCol.split('_')[0] + ':' + strType
            listDataColumnsClean.append(strColClean)

        # print('Collapsing microarray data to median values, this may take some time..')
        # arrayMedianProbeAbundances = np.zeros((numUniqueGenes, numMice), dtype=np.float)
        # for iGene in range(len(listUniqueMsENSEMBLGenes)):
        #     strENSGene = listUniqueMsENSEMBLGenes[iGene]
        #     listProbes = dfAnnot['PROBEID'][dfAnnot['ENSEMBL'] == strENSGene].tolist()
        #     arrayMedianProbeAbundances[iGene,:] = np.median(dfIn.reindex(listProbes).values.astype(np.float),
        #                                                     axis=0)

        dfData = pd.DataFrame(data=arrayMedianProbeAbundances,
                              index=listUniqueMsENSEMBLGenes,
                              columns=listDataColumnsClean)

        return dfData

    def mollaoglu_combined(flagResult=False,
                           strDataLoc='D:\\data\\KateSutherland',
                           strAnnotFileName='annotation.txt',
                           strTempFile='Mollaoglu_merged.pickle',
                           flagPerformExtraction=False):

        if not os.path.exists(os.path.join(strDataLoc, strTempFile)):
            flagPerformExtraction = True

        if flagPerformExtraction:

            dfAnnot = pd.read_table(os.path.join(strDataLoc, strAnnotFileName),
                                    sep='\t', quoting=2, header=0, index_col=None)

            dictENSMUSToHGNCHu = dict(zip(dfAnnot['ENSEMBL'].values.tolist(),
                                          dfAnnot['humanSymbol'].values.tolist()))

            dfRNAseq = PreProcess.mollaoglu_rnaseq()
            listRNAseqIndex = dfRNAseq.index.tolist()
            listRNAseqGenes = [strCol for strCol in listRNAseqIndex if 'ENSMUSG' in strCol]
            listRNASeqSamples = dfRNAseq.columns.tolist()

            dfMicroarray = PreProcess.mollaoglu_microarray()
            listMicroarrayIndex = dfMicroarray.index.tolist()
            listMicroarrayGenes = [strCol for strCol in listMicroarrayIndex if 'ENSMUSG' in strCol]
            listMicroarraySamples = dfMicroarray.columns.tolist()

            listSharedGenes = sorted(list(set(listRNAseqGenes).intersection(listMicroarrayGenes)))

            dfMerged = pd.concat([dfMicroarray.reindex(listSharedGenes), dfRNAseq.reindex(listSharedGenes)], axis=1)

            numMice = np.shape(dfMerged)[1]

            listDataHuGenes = [dictENSMUSToHGNCHu[strGene] for strGene in listSharedGenes]

            setDataHuGenes = set(listDataHuGenes)
            setDataHuGenes.remove(np.nan)

            # Load the NK cell transcriptomic signatures
            dictLocalNKSig = GeneSetScoring.CuratedList.nk_cells_in_tumour()
            listNKUpGenesInData = list(set(dictLocalNKSig['UpGenes']).intersection(setDataHuGenes))

            dictCombLM22TCellSig = GeneSetScoring.CuratedList.t_cells_in_tumour()
            listTCellUpGenesInData = list(set(dictCombLM22TCellSig['T cells']['UpGenes']).intersection(setDataHuGenes))

            dictOfDictLM22 = GeneSetScoring.ExtractList.cibersort_genes()
            listLM22Sets = list(dictOfDictLM22.keys())
            listLM22SetScoreHeaders = ['LM22: ' + strGeneSet + ' score' for strGeneSet in listLM22Sets]

            dfLM22Scores = pd.DataFrame(
                data=np.zeros((numMice, len(listLM22Sets))),
                columns=listLM22SetScoreHeaders,
                index=dfMerged.columns.tolist())

            print('Scoring mouse samples for specified signatures..')
            arrayNKScore = np.zeros(numMice, dtype=np.float)
            arrayTCellScore = np.zeros(numMice, dtype=np.float)

            for iMouse in range(numMice):
                print('{}'.format(iMouse + 1) + ' of ' + '{}'.format(numMice + 1))
                arrayNKScore[iMouse] = \
                    GeneSetScoring.FromInput.single_sample_rank_score(
                        listAllGenes=listDataHuGenes,
                        arrayTranscriptAbundance=dfMerged.iloc[:, iMouse].values.astype(np.float),
                        listUpGenesToScore=listNKUpGenesInData,
                        flagApplyNorm=True)
                arrayTCellScore[iMouse] = \
                    GeneSetScoring.FromInput.single_sample_rank_score(
                        listAllGenes=listDataHuGenes,
                        arrayTranscriptAbundance=dfMerged.iloc[:, iMouse].values.astype(np.float),
                        listUpGenesToScore=listTCellUpGenesInData,
                        flagApplyNorm=True)

                for strGeneSet in listLM22Sets:
                    listUpGenes = dictOfDictLM22[strGeneSet]['UpGenes']
                    dfLM22Scores['LM22: ' + strGeneSet + ' score'].iloc[
                        iMouse] = GeneSetScoring.FromInput.single_sample_rank_score(
                        listAllGenes=listDataHuGenes,
                        arrayTranscriptAbundance=dfMerged.iloc[:, iMouse].values.astype(np.float),
                        listUpGenesToScore=list(set(listUpGenes).intersection(setDataHuGenes)),
                        flagApplyNorm=True)

            dfOut = pd.concat([dfMerged.transpose().copy(), dfLM22Scores], axis=1)

            dfOut['NK score'] = arrayNKScore
            dfOut['T Cell score'] = arrayTCellScore

            dfOut.to_pickle(os.path.join(strDataLoc, strTempFile))

        else:

            dfOut = pd.read_pickle(os.path.join(strDataLoc, strTempFile))


        return dfOut

    def ms_hgnc_to_ensmusg(flagResult=False,
                           strDataLoc='D:\\data\\KateSutherland',
                           strAnnotFileName='annotation.txt'):

        dfAnnot = pd.read_table(os.path.join(strDataLoc, strAnnotFileName),
                                sep='\t', quoting=2, header=0, index_col=None)

        dictMsHGNCToENSMUSG = dict(zip(dfAnnot['mouseSymbol'].values.tolist(),
                                      dfAnnot['ENSEMBL'].values.tolist()))

        return dictMsHGNCToENSMUSG


    def ensmusg_to_hu_hgnc(flagResult=False,
                           strDataLoc='D:\\data\\KateSutherland',
                           strAnnotFileName='annotation.txt'):

        dfAnnot = pd.read_table(os.path.join(strDataLoc, strAnnotFileName),
                                sep='\t', quoting=2, header=0, index_col=None)

        dictENSMUSGToHuHGNC = dict(zip(dfAnnot['ENSEMBL'].values.tolist(),
                                      dfAnnot['humanSymbol'].values.tolist()))

        return dictENSMUSGToHuHGNC

    def rudin_classifications(
            flagResult=False,
            strFileName='41568_2019_164_MOESM1_ESM.xlsx'):

        dfIn = pd.read_excel(os.path.join(PathDir.strDataLoc, strFileName),
                             sheet_name='Sheet1', header=0)
        # for some reason the first row is NaN?
        dfIn.drop(labels=[0], axis=0, inplace=True)
        # dfIn.set_index(keus='Name', inplace=True)

        listSubtypes = dfIn['Subtype assignment'].unique().tolist()

        listIsCellLine = \
            dfIn['Name'][dfIn['Tumour or cell line']=='Cell line'].astype(str).tolist()

        listIsTumour = \
            dfIn['Name'][dfIn['Tumour or cell line']=='Tumour'].astype(str).tolist()

        dictSubtypes = {'Tumour':{},
                        'Cell line':{}}
        for strSubtype in listSubtypes:
            listInSubtype = dfIn['Name'][dfIn['Subtype assignment']==strSubtype].astype(str).tolist()

            dictSubtypes['Tumour'][strSubtype] = \
                sorted(list(set(listInSubtype).intersection(set(listIsTumour))))
            dictSubtypes['Cell line'][strSubtype] = \
                sorted(list(set(listInSubtype).intersection(set(listIsCellLine))))

        return dictSubtypes

    def george_subsets(flagResult=False):

        dfGeorgeSCLC = \
            PreProcess.george_2015_small_cell(flagPerformExtraction=False)
        listGeorgeSCLCPats = dfGeorgeSCLC.index.tolist()

        dictSubtypes = PreProcess.rudin_classifications()
        dictSubsets = dictSubtypes['Tumour']

        # NB: Patient sample classifications are taken directly from Rudin et al (2019), Nat Rev Cancer; Fig. S1
        listASCL1PosTumours = list(set(dictSubsets['SCLC-A']).intersection(set(listGeorgeSCLCPats)))
        listNEUROD1PosTumours = list(set(dictSubsets['SCLC-N']).intersection(set(listGeorgeSCLCPats)))
        listPOU2F3PosTumours = list(set(dictSubsets['SCLC-P']).intersection(set(listGeorgeSCLCPats)))
        listYAP1PosTumours = list(set(dictSubsets['SCLC-Y']).intersection(set(listGeorgeSCLCPats)))

        listToClassify = sorted(list(set(listGeorgeSCLCPats).difference(
            set(listASCL1PosTumours + listNEUROD1PosTumours + listPOU2F3PosTumours + listYAP1PosTumours))))

        for strSample in listToClassify:

            flagIsASCL1 = False
            flagIsNEUROD1 = False
            flagIsPOU2F3 = False
            flagIsYAP1 = False

            numASCL1Abund = np.log2(dfGeorgeSCLC.loc[strSample, 'ASCL1'].astype(np.float) + 1)
            numNEUROD1Abund = np.log2(dfGeorgeSCLC.loc[strSample, 'NEUROD1'].astype(np.float) + 1)
            numPOU2F3Abund = np.log2(dfGeorgeSCLC.loc[strSample, 'POU2F3'].astype(np.float) + 1)
            numYAP1Abund = np.log2(dfGeorgeSCLC.loc[strSample, 'YAP1'].astype(np.float) + 1)

            if numPOU2F3Abund > PreProcess.dictTFConstThresh['POU2F3']:
                flagIsPOU2F3 = True

            if numYAP1Abund > PreProcess.dictTFConstThresh['YAP1']:
                flagIsYAP1 = True

            if np.bitwise_and(numASCL1Abund > PreProcess.dictTFRelBaseThresh['ASCL1'],
                              (numASCL1Abund - PreProcess.dictTFRelDiffThresh['ASCL1']) > numNEUROD1Abund):
                flagIsASCL1 = True

            if np.bitwise_and(numNEUROD1Abund > PreProcess.dictTFRelBaseThresh['NEUROD1'],
                              (numNEUROD1Abund - PreProcess.dictTFRelDiffThresh['NEUROD1']) > numASCL1Abund):
                flagIsNEUROD1 = True

            numSummedFlags = np.sum([flagIsPOU2F3, flagIsYAP1, flagIsASCL1, flagIsNEUROD1])
            if numSummedFlags == 1:
                if flagIsPOU2F3:
                    listPOU2F3PosTumours.append(strSample)
                elif flagIsYAP1:
                    listYAP1PosTumours.append(strSample)
                elif flagIsASCL1:
                    listASCL1PosTumours.append(strSample)
                elif flagIsNEUROD1:
                    listNEUROD1PosTumours.append(strSample)

            elif numSummedFlags == 0:
                print(f'Warning: George::{strSample} satisfies no classification criteria and will remain unclassified')

            elif numSummedFlags == 2:
                print(f'Warning: George::{strSample} satisfies 2 classification criteria')
                if flagIsYAP1 & flagIsASCL1:
                    print(f'\t\t{strSample} matches both YAP1 & ASCL1 criteria --> excluding')
                elif flagIsPOU2F3 & flagIsASCL1:
                    print(f'\t\t{strSample} matches both POU2F3 & ASCL1 criteria --> excluding')
            else:
                print(f'Warning: George::{strSample} satisfies more than 2 classification criteria')

        listUnClassified = sorted(list(set(listGeorgeSCLCPats).difference(
            set(listASCL1PosTumours + listNEUROD1PosTumours + listPOU2F3PosTumours + listYAP1PosTumours))))

        # Unfortunately there was a version difference in the figures included with the paper
        #  submission that led to the exclusion of 4 samples which lacked metadata and/or definitive
        #  classifications
        listUnClassified.remove('S01297')
        listASCL1PosTumours.remove('S01563')
        listASCL1PosTumours.remove('S01248')
        listPOU2F3PosTumours.remove('S01556')

        return {'SCLC-A': listASCL1PosTumours,
                'SCLC-N': listNEUROD1PosTumours,
                'SCLC-P': listPOU2F3PosTumours,
                'SCLC-Y': listYAP1PosTumours,
                'Unclassified': listUnClassified}

    def george_stage(flagResult=False):

        # Load the metadata
        dfMeta = Load.george_2015_metadata()

        # Create a dictionary whhich maps index (sample ID) to clinical stage
        dictStage = dict(zip(dfMeta.index.tolist(),
                             dfMeta['stage_UICC'].values.tolist()))

        return dictStage

    def mouse_subsets(flagResult=False):

        strGEODataLoc = 'D:\\db\\geo\\GSE89660_-_Ms_SCLC'
        dfMsRNASeqAnnot = pd.read_table(os.path.join(strGEODataLoc, 'GSE89660_series_matrix.txt.gz'),
                                        skiprows=31, sep='\t', header=0, index_col=0)

        dictRNASeqNameToGenotype = dict(zip(dfMsRNASeqAnnot.columns.tolist(),
                                            dfMsRNASeqAnnot.iloc[9, :].values.tolist()))
        listRNASeqRPM = [strSmp for strSmp in dfMsRNASeqAnnot.columns.tolist()
                          if dictRNASeqNameToGenotype[strSmp] == 'genotype: RPM']
        listRNASeqRPR2 = [strSmp for strSmp in dfMsRNASeqAnnot.columns.tolist()
                          if dictRNASeqNameToGenotype[strSmp] == 'genotype: RPR2']

        dfMsData = PreProcess.mollaoglu_combined(flagPerformExtraction=False)
        listMice = dfMsData.index.tolist()

        listNormalMsLabels = [strSample for strSample in listMice if 'Normal' in strSample]
        listRPR2MsLabels = [strSample for strSample in listMice if 'TKO' in strSample] + listRNASeqRPR2
        listRPRMsLabels = [strSample for strSample in listMice if 'DKO' in strSample]
        listRPMLabels = listRNASeqRPM

        return {'RPM':listRPMLabels,
                'RPR':listRPRMsLabels,
                'RPR2':listRPR2MsLabels,
                'Normal':listNormalMsLabels}

class Plot:

    dictPlotColors = {'SCLC-A': '#0077b9',
                      'SCLC-N': '#6dbe83',
                      'RPR': '#0077b9',
                      'RPR2': '#0077b9',
                      'RPM': '#6dbe83',
                      'SCLC-P': '#004282',
                      'SCLC-Y': '#f26a6b',
                      'Unclassified': '#a9a3a1'}

    dictSubsetScatterMarker = {
        'SCLC-A': 'o',
        'SCLC-N': 'o',
        'SCLC-P': 's',
        'SCLC-Y': 'o',
        'Unclassified': 'o'}

    arrayColorMap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'SCLC-subsets',
        [(0 / 255, 119 / 255, 185 / 255, 1.0),      # SCLC-A
         (109 / 255, 190 / 255, 131 / 255, 1.0),    # SCLC-N
         (0, 30 / 255, 66 / 255, 1.0),              # SCLC-P
         (242 / 255, 106 / 255, 107 / 255, 1.0)],   # SCLC-Y
        N=4)

    arrayMouseColorMap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'SCLC-subsets',
        [(0 / 255, 140 / 255, 165 / 255, 1.0),      # RPR (~= SCLC-A)
         (30 / 255, 100 / 255, 205 / 255, 1.0),      # RPR2 (~= SCLC-A)
         (109 / 255, 190 / 255, 131 / 255, 1.0)],   # RPM (~= SCLC-N)
        N=3)

    listMsDataPlotOrder = ['RPR', 'RPR2', 'RPM']
    listHuDataPlotOrder = ['SCLC-A', 'SCLC-N', 'SCLC-P', 'SCLC-Y']

    numFontSize = 6

    listStageOrder = ['Unclassified', 'I', 'Ia', 'Ib', 'II', 'IIa', 'IIb', 'III', 'IIIa', 'IIIb', 'IV']
    listStageOrderLower = [strStage.lower() for strStage in listStageOrder]

    def figure_one(flagResult=False):

        listScoresForBoxPlots = ['NK score', 'T Cell score']
        # 'pSTAT1 score', 'DC score',

        listMarkerTFs = ['ASCL1', 'NEUROD1', 'POU2F3', 'YAP1']

        listOfHMGeneGroupLabels = ['MHC & antigen\npresentation',
                                   'Immune feedback.\n& regulation',
                                   'NK cell\ntargeting',
                                   'NK cell\nmarkers']
        listOfListsHMGenes = [['B2M', 'HLA-A',
                               'HLA-B', 'HLA-C', 'HLA-DMA', 'HLA-DMB',
                               'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1',
                               'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB2',
                               'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5',
                               'HLA-E', 'HLA-F', 'HLA-G'],
                              ['LAG3', 'FOXP3', 'IL6', 'CD86', 'CD96',
                               'TGFB1', 'IL10', 'CD80', 'ICOS', 'TIGIT',
                               'PDCD1', 'KLRK1'],
                              ['ULBP1', 'ULBP2', 'ULBP3', 'RAET1E', 'RAET1G', 'RAET1L',
                               'MICA', 'MICB'],
                              ['KLRC1', 'KLRC2', 'KLRC4', 'KLRD1', 'KLRF1',
                               'GZMB', 'XCL1', 'XCL2',
                               'NCR1', 'NCR3']]

        listHeatMapGenes = [strGene for listGenes in listOfListsHMGenes
                            for strGene in listGenes]

        listGenesForSwarm = ['CD274', 'TAP1', 'TAP2']
        arraySwarmGeneYTicks = np.arange(start=0, stop=9.1, step=3.0)

        numFigWidth = 5
        numFigHeight = 4

        numLeftLimit = 0.10
        numTopLimit = 0.96
        numBotLimit = 0.10
        numRightLimit = 0.90
        numSidePlotWidth = 0.15
        numHeatMapLeft = numLeftLimit+numSidePlotWidth+0.15
        numHeatMapWidth = numRightLimit-numHeatMapLeft

        arrayGridSpecSidePlots = matplotlib.gridspec.GridSpec(
            nrows=3, ncols=1, left=numLeftLimit, right=numLeftLimit+numSidePlotWidth,
            bottom=numBotLimit, top=numTopLimit,
            hspace=0.65)

        arraySubtypeHMPos = np.array([numHeatMapLeft, 0.94, numHeatMapWidth, 0.02],
                                           dtype=np.float)

        arrayMarkerTFHMPos = np.array([numHeatMapLeft, 0.87, numHeatMapWidth, 0.06],
                                           dtype=np.float)

        arrayStageHMPos = np.array([numHeatMapLeft, 0.84, numHeatMapWidth, 0.02],
                                           dtype=np.float)
        arrayStageHMCBarPos = np.array([numRightLimit+0.02, 0.79, 0.015, 0.14],
                                        dtype=np.float)

        arrayGeneAbundHMPos = np.array([numHeatMapLeft, numBotLimit, numHeatMapWidth, 0.73],
                                        dtype=np.float)
        arrayGeneAbundHMCBarPos = np.array([numRightLimit+0.02, 0.15, 0.015, 0.20],
                                        dtype=np.float)

        dfGeorgeSCLC = PreProcess.george_2015_small_cell(flagPerformExtraction=False)
        listPatients = dfGeorgeSCLC.index.tolist()

        dictSubsets = PreProcess.george_subsets()
        listSubsets = Plot.listHuDataPlotOrder

        dictStage = PreProcess.george_stage()

        listClassifiedPatients = \
            dictSubsets['SCLC-A'] + dictSubsets['SCLC-N'] + \
            dictSubsets['SCLC-P'] + dictSubsets['SCLC-Y']

        numPatients = len(listClassifiedPatients)

        listDataForSwarms = []
        for iGene in range(len(listGenesForSwarm)):
            arrayGeneAbund = dfGeorgeSCLC[listGenesForSwarm[iGene]].reindex(listClassifiedPatients).values.astype(np.float)
            listDataForSwarms.append(np.log2(arrayGeneAbund+1))

        arrayPatientClassHM = np.zeros((1,numPatients), dtype=np.float)
        listPatientsGrouped = []
        numPatientsClassified = 0
        listOfListsGrouped = []
        for iSubset in range(len(listSubsets)):
            strSubset = listSubsets[iSubset]
            listPatientsInGroup = list(set(dictSubsets[strSubset]).intersection(set(listClassifiedPatients)))
            numPatientsInSubset = len(listPatientsInGroup)
            arrayPatientClassHM[0,numPatientsClassified:numPatientsClassified+numPatientsInSubset] = np.float(iSubset)
            numPatientsClassified = numPatientsClassified + numPatientsInSubset
            listPatientsGrouped = listPatientsGrouped + listPatientsInGroup
            listOfListsGrouped.append(listPatientsInGroup)

        arrayTumorStageHM = np.zeros((1,numPatients), dtype=np.float)
        for iPatient in range(len(listPatientsGrouped)):
            strPatient = listPatientsGrouped[iPatient]
            if strPatient in dictStage.keys():
                strStage = dictStage[strPatient]
                if not (strStage == strStage):
                    strStage = 'Unclassified'
            else:
                strStage = 'Unclassified'

            numStage = Plot.listStageOrderLower.index(strStage.lower())
            arrayTumorStageHM[0,iPatient] = numStage

        dfStage = pd.DataFrame(data=arrayTumorStageHM,
                               columns=listPatientsGrouped,
                               index=['Stage'])

        arrayGeneAbundSorted = \
            np.log2(dfGeorgeSCLC[listHeatMapGenes].reindex(listPatientsGrouped).transpose().values.astype(np.float)+1)
        arrayGeneAbundZScoreNorm = np.zeros(np.shape(arrayGeneAbundSorted), dtype=np.float)
        for iGene in range(len(listHeatMapGenes)):
            numAbundMean = np.mean(arrayGeneAbundSorted[iGene,:])
            numAbundStDev = np.std(arrayGeneAbundSorted[iGene,:])
            arrayGeneAbundZScoreNorm[iGene,:] = (arrayGeneAbundSorted[iGene,:]-numAbundMean)/numAbundStDev

        dfGeneAbundZNorm = pd.DataFrame(data=arrayGeneAbundZScoreNorm,
                                        index=listHeatMapGenes,
                                        columns=listPatientsGrouped)

        arrayMarkerTFSorted = \
            np.log2(dfGeorgeSCLC[listMarkerTFs].reindex(listPatientsGrouped).transpose().values.astype(np.float)+1)
        arrayMarkerTFZScoreNorm = np.zeros(np.shape(arrayMarkerTFSorted), dtype=np.float)
        for iGene in range(len(listMarkerTFs)):
            numAbundMean = np.mean(arrayMarkerTFSorted[iGene,:])
            numAbundStDev = np.std(arrayMarkerTFSorted[iGene,:])
            arrayMarkerTFZScoreNorm[iGene,:] = (arrayMarkerTFSorted[iGene,:]-numAbundMean)/numAbundStDev

        dfMarkerTFZNorm = pd.DataFrame(data=arrayMarkerTFZScoreNorm,
                                        index=listMarkerTFs,
                                        columns=listPatientsGrouped)

        listPatientPlottingOrder = []
        for iSubset in range(len(listSubsets)):
            listPatientsInSubset = listOfListsGrouped[iSubset]

            arrayStageOrder = np.argsort(dfStage[listPatientsInSubset].loc['Stage'].values.astype(np.float))
            listPatientClust = [listPatientsInSubset[i] for i in arrayStageOrder]

            listPatientPlottingOrder = listPatientPlottingOrder + listPatientClust

        handFig = plt.figure()
        handFig.set_size_inches(w=numFigWidth,h=numFigHeight)

        handAx = plt.subplot(arrayGridSpecSidePlots[0])

        sns.swarmplot(data=listDataForSwarms,
                      size=1.5,
                      color='0.5')

        handAx.set_xticks([])
        handAx.spines['right'].set_visible(False)
        handAx.spines['top'].set_visible(False)

        arrayYLim = handAx.get_ylim()
        numTextYPos = arrayYLim[0]-0.03*np.ptp(arrayYLim)

        for iGene in range(len(listGenesForSwarm)):
            handAx.text(iGene, numTextYPos,
                        listGenesForSwarm[iGene],
                        ha='center', va='top', rotation=90,
                        fontstyle='italic', fontsize=Plot.numFontSize*0.8)

        handAx.set_ylim([arraySwarmGeneYTicks[0]-0.2, arraySwarmGeneYTicks[-1]])
        handAx.set_yticks(arraySwarmGeneYTicks)
        for handTick in handAx.yaxis.get_major_ticks():
            handTick.label.set_fontsize(Plot.numFontSize)

        handAx.set_ylabel('log$_{2}$|FPKM+1|',
                          fontsize=Plot.numFontSize)

        for iScore in range(len(listScoresForBoxPlots)):
            strScore = listScoresForBoxPlots[iScore]

            listScoresToPlot = []
            for strSubset in listSubsets:
                listPatientsInSubset = list(set(dictSubsets[strSubset]).intersection(set(listPatients)))
                listScoresToPlot.append(
                    dfGeorgeSCLC[strScore].reindex(listPatientsInSubset).values.astype(np.float))

            handAx = plt.subplot(arrayGridSpecSidePlots[iScore+1])

            handBoxPlot = handAx.boxplot(listScoresToPlot, widths=0.2, notch=0, sym='+', vert=1, whis=1.5)
            plt.setp(handBoxPlot['boxes'], color='black', linewidth=0.5)
            plt.setp(handBoxPlot['whiskers'], color='black', linewidth=0.5)
            plt.setp(handBoxPlot['caps'], color='black', linewidth=0.5)
            plt.setp(handBoxPlot['medians'], color='orange', linewidth=0.5)
            plt.setp(handBoxPlot['fliers'], marker='')
            for iSubset in range(len(listSubsets)):
                strSubset = listSubsets[iSubset]
                plt.scatter(np.ones(len(listScoresToPlot[iSubset]), dtype=np.float) * (0.7 + np.float(iSubset)),
                            listScoresToPlot[iSubset],
                            color=Plot.dictPlotColors[strSubset],
                            s=3,
                            lw=0.0,
                            alpha=0.5)

            handAx.set_xticklabels([])
            arrayYLim = handAx.get_ylim()
            for i in range(len(listSubsets)):
                handAx.text(i+1.0, np.min(arrayYLim)-0.08*np.ptp(arrayYLim),
                            listSubsets[i].split('SCLC-')[1], fontsize=Plot.numFontSize, rotation=0,
                            ha='center', va='top')

            for handTick in handAx.yaxis.get_major_ticks():
                handTick.label.set_fontsize(Plot.numFontSize)
            handAx.set_title(strScore, fontsize=Plot.numFontSize)

            handAx.spines['right'].set_visible(False)
            handAx.spines['top'].set_visible(False)


        handAx = handFig.add_axes(arraySubtypeHMPos)
        handAx.matshow(arrayPatientClassHM,
                       cmap=Plot.arrayColorMap,
                       aspect='auto')
        handAx.set_xticks([])
        handAx.set_yticks([])

        handAx.text(-2, 0.0, 'Subtype', ha='right', va='center', fontsize=Plot.numFontSize)
        arrayYLim = handAx.get_ylim()
        numBase = 0
        for iSubtype in range(len(listSubsets)):
            numInSubset = len(listOfListsGrouped[iSubtype])
            handAx.text(numBase+(numInSubset/2), arrayYLim[1]+0.03*np.ptp(arrayYLim),
                        listSubsets[iSubtype].split('SCLC-')[1], ha='center', va='bottom',
                        fontsize=Plot.numFontSize)
            numBase = numBase + numInSubset

        handAx = handFig.add_axes(arrayMarkerTFHMPos)
        handAbundHM = handAx.matshow(dfMarkerTFZNorm[listPatientPlottingOrder].values.astype(np.float),
                       cmap=plt.cm.PRGn,
                       vmin=-3.0,
                       vmax=3.0,
                       aspect='auto')
        handAx.set_xticks([])
        handAx.set_yticks([])
        for iGene in range(len(listMarkerTFs)):
            handAx.text(-1.5, iGene, listMarkerTFs[iGene],
                        ha='right', va='center', fontsize=Plot.numFontSize*0.75,
                        style='italic')

        structAxPos = handAx.get_position()

        numGenesInSubset = len(listOfListsHMGenes[iSubset])
        numTextYPos = structAxPos.y0 + 0.5*structAxPos.height
        handFig.text(numHeatMapLeft-0.10, numTextYPos, 'Marker\nTFs',
                     rotation=90, fontsize=Plot.numFontSize*0.75, weight='bold', ha='center', va='center')

        handAx = handFig.add_axes(arrayStageHMPos)
        handStageHM = handAx.matshow(dfStage[listPatientPlottingOrder],
                       cmap=plt.get_cmap('plasma', len(Plot.listStageOrder)),
                       aspect='auto')
        handAx.set_xticks([])
        handAx.set_yticks([])

        handAx.text(-2, 0.0, 'Stage', ha='right', va='center', fontsize=Plot.numFontSize)

        handAxStageCBar = handFig.add_axes(arrayStageHMCBarPos)
        handColorBar = plt.colorbar(handStageHM, cax=handAxStageCBar, ticks=[])
        structAxPos = handAxStageCBar.get_position()
        numRowHeight = structAxPos.height/len(Plot.listStageOrder)
        for iStage in range(len(Plot.listStageOrder)):
            handFig.text(structAxPos.x0+1.2*structAxPos.width,
                         structAxPos.y0+(iStage+0.5)*numRowHeight,
                         Plot.listStageOrder[iStage],
                         ha='left', va='center',
                         fontsize=Plot.numFontSize*0.5)

        # handColorBar.ax.tick_params(labelsize=Plot.numFontSize)

        handAx = handFig.add_axes(arrayGeneAbundHMPos)
        handAbundHM = handAx.matshow(dfGeneAbundZNorm[listPatientPlottingOrder].values.astype(np.float),
                       cmap=plt.cm.PRGn,
                       vmin=-3.0,
                       vmax=3.0,
                       aspect='auto')
        handAx.set_xticks([])
        handAx.set_yticks([])
        for iGene in range(len(listHeatMapGenes)):
            handAx.text(-1.5, iGene, listHeatMapGenes[iGene],
                        ha='right', va='center', fontsize=Plot.numFontSize*0.75,
                        style='italic')

        structAxPos = handAx.get_position()
        numGroupedGenesLabelled = 0
        numHeightPerRow = structAxPos.height / len(listHeatMapGenes)
        for iSubset in range(len(listOfHMGeneGroupLabels)):
            strLabel = listOfHMGeneGroupLabels[iSubset]
            numGenesInSubset = len(listOfListsHMGenes[iSubset])
            numTextYPos = (structAxPos.y0 + structAxPos.height) - \
                          (numGroupedGenesLabelled+(numGenesInSubset/2))*numHeightPerRow
            handFig.text(numHeatMapLeft-0.10, numTextYPos, strLabel,
                         rotation=90, fontsize=Plot.numFontSize*0.75, weight='bold', ha='center', va='center')
            if iSubset < len(listOfHMGeneGroupLabels):
                handAx.axhline(y=numGroupedGenesLabelled +numGenesInSubset-0.5, xmin=0.0, xmax=1.0, color='w', lw=1.0)
                handAx.axhline(y=numGroupedGenesLabelled +numGenesInSubset-0.5, xmin=0.0, xmax=1.0, color='k', lw=0.5)
            numGroupedGenesLabelled = numGroupedGenesLabelled + numGenesInSubset

        handAxCBar = handFig.add_axes(arrayGeneAbundHMCBarPos)
        handColorBar = plt.colorbar(handAbundHM, cax=handAxCBar)
        handColorBar.ax.tick_params(labelsize=Plot.numFontSize*0.7)
        handAxCBar.set_title('$z$-score', fontsize=Plot.numFontSize*0.7)

        handFig.savefig(os.path.join(PathDir.strOutputLoc, 'Figure1_CDEF.png'),
                        ext='png', dpi=300)
        handFig.savefig(os.path.join(PathDir.strOutputLoc, 'Figure1_CDEF.pdf'),
                        ext='pdf', dpi=300)
        plt.close(handFig)

        return flagResult

    def figure_two(flagResult=False):

        listScoresForBoxPlots = ['NK score', 'T Cell score']

        listOfHMGeneGroupLabels = ['Marker\nTFs',
                                   'MHC & antigen\npresentation',
                                   'Immune feedback.\n& regulation',
                                   'NK cell\ntargeting',
                                   'NK cell\nmarkers']
        listOfListsHMGenes = [['Ascl1', 'Neurod1', 'Myc'],
                              ['B2m', 'H2-T23', 'H2-M3',
                               'H2-Eb1', 'H2-Eb2', 'H2-Aa'],
                              # NOT PRESENT: 'H2-D1','H2-Q10','H2-Q2','H2-Q7'
                              ['Lag3', 'Cd274', 'Foxp3', 'Il6', 'Cd86', 'Cd96',
                               'Tgfb1', 'Il10', 'Cd80', 'Icos',
                               'Pdcd1', 'Klrk1'],
                              ['Ulbp1', 'Raet1d'],
                              ['Ctsw', 'Fasl', 'Gzmb', 'Klrb1', 'Klrc1', 'Klrc2', 'Klrd1',
                               'Ncr1', 'Nkg7', 'Prf1', 'Xcl1']
                              # NOT PRESENT: 'Klrc4', 'Klrf1', 'Xcl2', 'Ncr3'
                              ]

        listHeatMapGenes = [strGene for listGenes in listOfListsHMGenes
                            for strGene in listGenes]

        numFigWidth = 5
        numFigHeight = 6

        numLeftLimit = 0.14
        numRightLimit = 0.87

        arrayDataSourceHMPos = np.array([numLeftLimit, 0.77, (numRightLimit-numLeftLimit), 0.02],
                                           dtype=np.float)

        arraySubsetHMPos = np.array([numLeftLimit, 0.74, (numRightLimit-numLeftLimit), 0.02],
                                           dtype=np.float)

        arrayGeneAbundHMPos = np.array([numLeftLimit, 0.22, (numRightLimit-numLeftLimit), 0.51],
                                        dtype=np.float)
        arrayGeneAbundHMCBarPos = np.array([numRightLimit+0.05, 0.32+(0.38/4), 0.02, 0.38/2],
                                        dtype=np.float)

        arrayGridSpecScores = matplotlib.gridspec.GridSpec(
            nrows=1, ncols=len(listScoresForBoxPlots), left=numLeftLimit, right=numRightLimit,
            bottom=0.10, top=0.21,
            wspace=0.90)

        dfMsData = PreProcess.mollaoglu_combined(flagPerformExtraction=False)
        listMice = dfMsData.index.tolist()

        dictSubsets = PreProcess.mouse_subsets()

        dictMiceToDataType = dict(zip(listMice, [None]*len(listMice)))
        for strMouse in listMice:
            if np.bitwise_and(strMouse[0:len('GSM')] == 'GSM',
                              ':' in strMouse):
                dictMiceToDataType[strMouse] = 'Microarray'
            elif np.bitwise_or(strMouse[0:len('PB')] == 'PB',
                               strMouse[0:len('S')] == 'S'):
                dictMiceToDataType[strMouse] = 'RNA-seq'

        listMousePlottingOrder = []
        for strSubset in Plot.listMsDataPlotOrder:
            listMiceInSubset = list(set(dictSubsets[strSubset]).intersection(set(listMice)))
            listSubsetMicroarraySamples = [strSample for strSample in listMiceInSubset
                                           if dictMiceToDataType[strSample]=='Microarray']
            listSubsetRNAseqSamples = [strSample for strSample in listMiceInSubset
                                       if dictMiceToDataType[strSample]=='RNA-seq']

            listMousePlottingOrder = listMousePlottingOrder + \
                                     listSubsetMicroarraySamples + listSubsetRNAseqSamples

        numMiceToPlot = len(listMousePlottingOrder)

        dictMsHGNCToENSMUSG = PreProcess.ms_hgnc_to_ensmusg()

        # fix some incorrect identifiers tot ry and improve coverage
        dictMsHGNCToENSMUSG['Dlk1'] = 'ENSMUSG00000040856'
        dictMsHGNCToENSMUSG['Ulbp1'] = 'ENSMUSG00000079685'
        dictMsHGNCToENSMUSG['Gzmb'] = 'ENSMUSG00000015437'
        dictMsHGNCToENSMUSG['Klrb1'] = 'ENSMUSG00000030361'
        dictMsENSMUSGToHGNC = dict(zip(dictMsHGNCToENSMUSG.values(), dictMsHGNCToENSMUSG.keys()))
        dictMsENSMUSGToHGNC['ENSMUSG00000032946'] = 'Rasgrp2'
        dictMsENSMUSGToHGNC['ENSMUSG00000071068'] = 'Treml2'

        listDataCols = dfMsData.columns.tolist()

        listGenes = [strCol for strCol in listDataCols
                     if strCol[0:len('ENSMUSG')] == 'ENSMUSG']

        for strENSMUSG in listGenes:
            if strENSMUSG not in dictMsHGNCToENSMUSG.keys():
                dictMsHGNCToENSMUSG['ENSMUSG'] = 'failed_map'

        dfPercentiles = pd.DataFrame(data=np.zeros((numMiceToPlot, len(listHeatMapGenes)), dtype=np.float),
                                     columns=listHeatMapGenes,
                                     index=listMousePlottingOrder)

        listHeatMapGenesENSMUSG = [dictMsHGNCToENSMUSG[strGene] for strGene in listHeatMapGenes]

        for iMouse in range(numMiceToPlot):
            strMouse = listMousePlottingOrder[iMouse]
            arrayAllAbundData = dfMsData.loc[strMouse].values.astype(np.float)
            arraySelGeneAbund = dfMsData[listHeatMapGenesENSMUSG].loc[strMouse].values.astype(np.float)
            for iGene in range(len(listHeatMapGenesENSMUSG)):
                dfPercentiles.iloc[iMouse,iGene] = scipy.stats.percentileofscore(
                    arrayAllAbundData,
                    arraySelGeneAbund[iGene])

        dfPercentileZScore = pd.DataFrame(data=np.zeros((numMiceToPlot, len(listHeatMapGenes)), dtype=np.float),
                                     columns=listHeatMapGenes,
                                     index=listMousePlottingOrder)
        for iGene in range(len(listHeatMapGenes)):
            arrayGenePercentileAbund = dfPercentiles.iloc[:,iGene].values.astype(np.float)
            numPercAbundMean = np.mean(arrayGenePercentileAbund)
            numPercAbundStDev = np.std(arrayGenePercentileAbund)
            dfPercentileZScore.iloc[:,iGene] = (arrayGenePercentileAbund - numPercAbundMean)/numPercAbundStDev

        listMiceDataSources = sorted(list(set(dictMiceToDataType.values())))

        arrayDataSource = np.zeros((1,len(listMousePlottingOrder)),
                                   dtype=np.float)
        for iMouse in range(len(listMousePlottingOrder)):
            strMouse = listMousePlottingOrder[iMouse]
            strSource = dictMiceToDataType[strMouse]
            arrayDataSource[0,iMouse] = listMiceDataSources.index(strSource)+5

        handFig = plt.figure()
        handFig.set_size_inches(w=numFigWidth,h=numFigHeight)

        arrayMouseClassHM = np.zeros((1,numMiceToPlot), dtype=np.float)
        for iMouse in range(len(listMousePlottingOrder)):
            strMouse = listMousePlottingOrder[iMouse]
            for iSubset in range(len(Plot.listMsDataPlotOrder)):
                if strMouse in dictSubsets[Plot.listMsDataPlotOrder[iSubset]]:
                    arrayMouseClassHM[0,iMouse] = iSubset

        arrayColorNorm = matplotlib.colors.Normalize(vmin=0, vmax=9)

        arrayColorMap = matplotlib.cm.ScalarMappable(norm=arrayColorNorm, cmap=plt.cm.tab10)

        handAx = handFig.add_axes(arrayDataSourceHMPos)
        plt.scatter(-3, -3, s=2, marker='o', color=arrayColorMap.to_rgba(5), label='Microarray')
        plt.scatter(-3, -3, s=2, marker='o', color=arrayColorMap.to_rgba(6), label='RNA-seq')
        handAx.matshow(arrayDataSource,
                       cmap=plt.cm.tab10,
                       vmin=0, vmax=9,
                       aspect='auto')
        handAx.set_xlim([-0.5, np.shape(arrayDataSource)[1]-0.5])
        handAx.set_ylim([-0.5, np.shape(arrayDataSource)[0] - 0.5])
        handAx.set_xticks([])
        handAx.set_yticks([])

        plt.legend(loc='center left',
                   bbox_to_anchor=(1.01, 0.5),
                   fontsize=Plot.numFontSize * 0.5,
                   scatterpoints=1,
                   framealpha=1,
                   ncol=1)

        handAx.text(-1, 0.0, 'Data Source', ha='right', va='center', fontsize=Plot.numFontSize)

        handAx = handFig.add_axes(arraySubsetHMPos)
        plt.scatter(-3, -3, s=2, marker='o', color=Plot.arrayMouseColorMap(0), label=Plot.listMsDataPlotOrder[0])
        plt.scatter(-3, -3, s=2, marker='o', color=Plot.arrayMouseColorMap(1), label=Plot.listMsDataPlotOrder[1])
        plt.scatter(-3, -3, s=2, marker='o', color=Plot.arrayMouseColorMap(2), label=Plot.listMsDataPlotOrder[2])
        handAx.matshow(arrayMouseClassHM,
                       cmap=Plot.arrayMouseColorMap,
                       aspect='auto')
        handAx.set_xlim([-0.5, np.shape(arrayMouseClassHM)[1]-0.5])
        handAx.set_ylim([-0.5, np.shape(arrayMouseClassHM)[0] - 0.5])
        handAx.set_xticks([])
        handAx.set_yticks([])
        plt.legend(loc='upper left',
                   bbox_to_anchor=(1.01, 1.01),
                   fontsize=Plot.numFontSize * 0.5,
                   scatterpoints=1,
                   framealpha=1,
                   ncol=1)

        handAx.text(-1, 0.0, 'Subset', ha='right', va='center', fontsize=Plot.numFontSize)

        handAx = handFig.add_axes(arrayGeneAbundHMPos)
        handAbundHM = handAx.matshow(dfPercentileZScore[listHeatMapGenes].reindex(listMousePlottingOrder).transpose().values.astype(np.float),
                       cmap=plt.cm.PRGn,
                       vmin=-3,
                       vmax=3,
                       aspect='auto')
        handAx.set_xticks([])
        handAx.set_yticks([])
        for iGene in range(len(listHeatMapGenes)):
            handAx.text(-1.0, iGene, listHeatMapGenes[iGene],
                        ha='right', va='center', fontsize=Plot.numFontSize*0.75,
                        style='italic')

        structAxPos = handAx.get_position()
        numGroupedGenesLabelled = 0
        numHeightPerRow = structAxPos.height / len(listHeatMapGenes)
        for iSubset in range(len(listOfHMGeneGroupLabels)):
            strLabel = listOfHMGeneGroupLabels[iSubset]
            numGenesInSubset = len(listOfListsHMGenes[iSubset])
            numTextYPos = (structAxPos.y0 + structAxPos.height) - \
                          (numGroupedGenesLabelled+(numGenesInSubset/2))*numHeightPerRow
            handFig.text(0.04, numTextYPos, strLabel,
                         rotation=90, fontsize=Plot.numFontSize*0.75, weight='bold', ha='center', va='center')
            if iSubset < len(listOfHMGeneGroupLabels):
                handAx.axhline(y=numGroupedGenesLabelled +numGenesInSubset-0.5, xmin=0.0, xmax=1.0, color='w', lw=1.0)
                handAx.axhline(y=numGroupedGenesLabelled +numGenesInSubset-0.5, xmin=0.0, xmax=1.0, color='k', lw=0.5)
            numGroupedGenesLabelled = numGroupedGenesLabelled + numGenesInSubset

        handAxCBar = handFig.add_axes(arrayGeneAbundHMCBarPos)
        handColorBar = plt.colorbar(handAbundHM, cax=handAxCBar)
        handColorBar.ax.tick_params(labelsize=Plot.numFontSize)
        handAxCBar.set_title('Percentile\nabundance\n$z$-score', fontsize=Plot.numFontSize)

        for iScore in range(len(listScoresForBoxPlots)):
            strScore = listScoresForBoxPlots[iScore]

            listScoresToPlot = []
            for strSubset in Plot.listMsDataPlotOrder:
                listMiceInSubset = list(set(dictSubsets[strSubset]).intersection(set(listMice)))
                listScoresToPlot.append(
                    dfMsData[strScore].reindex(listMiceInSubset).values.astype(np.float))

            handAx = plt.subplot(arrayGridSpecScores[iScore])

            handBoxPlot = handAx.boxplot(listScoresToPlot, widths=0.2, notch=0, sym='+', vert=1, whis=1.5)
            plt.setp(handBoxPlot['boxes'], color='black', linewidth=0.5)
            plt.setp(handBoxPlot['whiskers'], color='black', linewidth=0.5)
            plt.setp(handBoxPlot['caps'], color='black', linewidth=0.5)
            plt.setp(handBoxPlot['medians'], color='orange', linewidth=0.5)
            plt.setp(handBoxPlot['fliers'], marker='')
            for iSubset in range(len(Plot.listMsDataPlotOrder)):
                strSubset = Plot.listMsDataPlotOrder[iSubset]
                plt.scatter(np.ones(len(listScoresToPlot[iSubset]), dtype=np.float) * (0.7 + np.float(iSubset)),
                            listScoresToPlot[iSubset],
                            color=Plot.dictPlotColors[strSubset],
                            s=3,
                            lw=0.0,
                            alpha=0.5)

            handAx.set_xticklabels([])
            arrayYLim = handAx.get_ylim()
            for i in range(len(Plot.listMsDataPlotOrder)):
                handAx.text(i+1.0, np.min(arrayYLim)-0.08*np.ptp(arrayYLim),
                            Plot.listMsDataPlotOrder[i], fontsize=Plot.numFontSize, rotation=90,
                            ha='center', va='top')

            for handTick in handAx.yaxis.get_major_ticks():
                handTick.label.set_fontsize(Plot.numFontSize)
            handAx.set_ylabel(strScore, fontsize=Plot.numFontSize)

            handAx.spines['right'].set_visible(False)
            handAx.spines['top'].set_visible(False)

        handFig.savefig(os.path.join(PathDir.strOutputLoc, 'Figure2.png'),
                        ext='png', dpi=300)
        handFig.savefig(os.path.join(PathDir.strOutputLoc, 'Figure2.pdf'),
                        ext='pdf', dpi=300)
        plt.close(handFig)

        return flagResult

    def supp_fig_one_scatter(flagResult=False):

        numScatterXTicks = 4
        numScatterYTicks = 4

        listGenePlottingOrder = ['NEUROD1', 'YAP1', 'POU2F3', 'ASCL1']
        listTopHistGenes = listGenePlottingOrder[-2:]
        listSideHistGenes = listGenePlottingOrder[0:2]

        arrayTopHistFreqLim = [0, 33]
        arraySideHistFreqLim = [0, 23]

        numFigWidth = 5
        numFigHeight = 5

        numTopHistHeight = 0.10
        numSideHistWidth = numTopHistHeight*(numFigHeight/numFigWidth)
        numHistXSpacing = 0.03
        numHistYSpacing = numHistXSpacing*(numFigWidth/numFigHeight)

        numLeftLimit = 0.10
        numRightLimit = 0.99
        numTopLimit = 0.99

        numScatterRightLimit = numRightLimit - (numSideHistWidth+numHistXSpacing)
        numScatterPlotWidthSpacing = 0.15
        numScatterPlotHeightSpacing = numScatterPlotWidthSpacing*(numFigHeight/numFigWidth)

        numCols = len(listGenePlottingOrder)-1
        numRows = len(listGenePlottingOrder)-1
        arrayGridSpecScatter = matplotlib.gridspec.GridSpec(
            nrows=numRows, ncols=numCols, left=numLeftLimit, right=numScatterRightLimit,
            bottom=0.09, top=numTopLimit - (numTopHistHeight+numHistYSpacing),
            hspace=numScatterPlotHeightSpacing, wspace=numScatterPlotWidthSpacing)


        # arrayed gene pairs; x,y
        listClassifyingPairsToComp = []
        for yGene in range(numCols):
            for xGene in range(yGene + 1, len(listGenePlottingOrder)):
                listClassifyingPairsToComp.append(
                    [listGenePlottingOrder[xGene], listGenePlottingOrder[yGene]])

        dfGeorgeSCLC = PreProcess.george_2015_small_cell(flagPerformExtraction=False)
        listPatients = dfGeorgeSCLC.index.tolist()

        dictSubsets = PreProcess.george_subsets()
        listSubsets = sorted(list(dictSubsets.keys()))

        dictGenePlotLim = {}
        for strGene in listGenePlottingOrder:
            arrayGeneAbund = np.log2(dfGeorgeSCLC[strGene].values.astype(np.float) + 1)
            numMinLim = np.min(arrayGeneAbund) - 0.05*np.ptp(arrayGeneAbund)
            numMaxLim = np.max(arrayGeneAbund) + 0.05*np.ptp(arrayGeneAbund)
            dictGenePlotLim[strGene] = (numMinLim, numMaxLim)

        handFig = plt.figure()
        handFig.set_size_inches(w=numFigWidth,h=numFigHeight)

        iRow = 0
        iCol = 0
        for iPair in range(len(listClassifyingPairsToComp)):
            strGeneOne = listClassifyingPairsToComp[iPair][0]
            strGeneTwo = listClassifyingPairsToComp[iPair][1]

            handAx = plt.subplot(arrayGridSpecScatter[iRow, iCol])
            for iGroup in range(len(listSubsets)):
                strGroup = listSubsets[iGroup]
                listPatientsInGroup = list(set(dictSubsets[strGroup]).intersection(set(listPatients)))
                handAx.scatter(np.log2(dfGeorgeSCLC[strGeneOne].reindex(listPatientsInGroup) + 1),
                               np.log2(dfGeorgeSCLC[strGeneTwo].reindex(listPatientsInGroup) + 1),
                               s=5,
                               marker=Plot.dictSubsetScatterMarker[strGroup],
                               alpha=0.7, lw=0.0,
                               c=Plot.dictPlotColors[strGroup],
                               label=strGroup)

            handAx.set_xlim([dictGenePlotLim[strGeneOne][0], dictGenePlotLim[strGeneOne][1]])
            arrayXTickLoc = plt.MaxNLocator(numScatterXTicks)
            handAx.xaxis.set_major_locator(arrayXTickLoc)

            handAx.set_ylim([dictGenePlotLim[strGeneTwo][0], dictGenePlotLim[strGeneTwo][1]])
            arrayYTickLoc = plt.MaxNLocator(numScatterYTicks)
            handAx.yaxis.set_major_locator(arrayYTickLoc)

            handAx.spines['right'].set_visible(False)
            handAx.spines['top'].set_visible(False)

            for axis in ['bottom', 'left']:
                handAx.spines[axis].set_linewidth(0.5)
            handAx.tick_params(width=0.5)

            if strGeneOne in PreProcess.dictTFConstThresh.keys():
                handAx.axvline(x=PreProcess.dictTFConstThresh[strGeneOne],
                               ymin=0.0, ymax=1.0,
                               linestyle='--', linewidth=0.5, color='0.5', alpha=0.5)

            if strGeneTwo in PreProcess.dictTFConstThresh.keys():
                handAx.axhline(y=PreProcess.dictTFConstThresh[strGeneTwo],
                               xmin=0.0, xmax=1.0,
                               linestyle='--', linewidth=0.5, color='0.5', alpha=0.5)

            if np.bitwise_and(strGeneOne == 'ASCL1', strGeneTwo == 'NEUROD1'):
                arrayXLim = handAx.get_xlim()
                arrayYLim = handAx.get_ylim()
                numMinMax = np.min([arrayXLim[1], arrayYLim[1]])

                handAx.plot([PreProcess.dictTFRelBaseThresh['ASCL1'], PreProcess.dictTFRelBaseThresh['ASCL1']],
                            [arrayXLim[0], PreProcess.dictTFRelBaseThresh['ASCL1'] - PreProcess.dictTFRelDiffThresh['NEUROD1']],
                            linestyle='--', linewidth=0.5, color='0.5', alpha=0.5)
                handAx.plot([PreProcess.dictTFRelBaseThresh['ASCL1'], numMinMax+PreProcess.dictTFRelDiffThresh['NEUROD1']],
                            [PreProcess.dictTFRelBaseThresh['ASCL1'] - PreProcess.dictTFRelDiffThresh['NEUROD1'],
                             arrayYLim[1]],
                            linestyle='--', linewidth=0.5, color='0.5', alpha=0.5)

                handAx.plot([arrayXLim[0], PreProcess.dictTFRelBaseThresh['NEUROD1'] - PreProcess.dictTFRelDiffThresh['NEUROD1']],
                            [PreProcess.dictTFRelBaseThresh['NEUROD1'], PreProcess.dictTFRelBaseThresh['NEUROD1']],
                            linestyle='--', linewidth=0.5, color='0.5', alpha=0.5)
                handAx.plot([PreProcess.dictTFRelBaseThresh['NEUROD1'] - PreProcess.dictTFRelDiffThresh['NEUROD1'],
                             numMinMax - PreProcess.dictTFRelDiffThresh['NEUROD1']],
                            [PreProcess.dictTFRelBaseThresh['NEUROD1'], arrayYLim[1]],
                            linestyle='--', linewidth=0.5, color='0.5', alpha=0.5)


                handAx.plot([0, numMinMax], [0, numMinMax],
                            linestyle='--', linewidth=0.5, color='k', alpha=0.5)
                handAx.text(numMinMax*1.00, numMinMax*1.01,
                            '$x=y$',
                            ha='center', va='bottom',
                            fontsize=Plot.numFontSize*0.5)

            if iRow == iCol:
                handAx.set_xlabel('$'+strGeneOne+'$' + '\nlog$_{2}$(FPKM+1)',
                                  fontsize=Plot.numFontSize)
                handAx.set_ylabel('$'+strGeneTwo +'$' + '\nlog$_{2}$(FPKM+1)',
                                  fontsize=Plot.numFontSize)

                for handTick in handAx.xaxis.get_major_ticks():
                    handTick.label.set_fontsize(Plot.numFontSize)

                for handTick in handAx.yaxis.get_major_ticks():
                    handTick.label.set_fontsize(Plot.numFontSize)

            else:
                handAx.set_xticklabels([])
                handAx.set_yticklabels([])

            if np.bitwise_and(iRow == 1, iCol == 1):
                plt.legend(loc='upper left',
                           bbox_to_anchor=(-1.4, 0.7),
                           fontsize=Plot.numFontSize*1.25,
                           scatterpoints=1,
                           framealpha=1,
                           ncol=1)

                structAxPos = handAx.get_position()
                numHorizHistXStart = structAxPos.xmin
                numVertHistYStart = structAxPos.ymin

            iCol += 1
            if iCol >= numCols:
                iRow += 1
                iCol = iRow

        arrayGridSpecTopHist = matplotlib.gridspec.GridSpec(
            nrows=1, ncols=2, left=numHorizHistXStart, right=numScatterRightLimit,
            bottom=numTopLimit - numTopHistHeight, top=numTopLimit,
            hspace=numScatterPlotHeightSpacing, wspace=numScatterPlotWidthSpacing)

        arrayGridSpecSideHist = matplotlib.gridspec.GridSpec(
            nrows=2, ncols=1, left=numRightLimit-numSideHistWidth, right=numRightLimit,
            bottom=numVertHistYStart, top=numTopLimit - (numTopHistHeight+numHistYSpacing),
            hspace=numScatterPlotHeightSpacing, wspace=numScatterPlotWidthSpacing)

        for iGene in range(len(listTopHistGenes)):

            arrayHistFreqTicks = np.arange(start=arrayTopHistFreqLim[0], stop=arrayTopHistFreqLim[1], step=10)

            strGene = listTopHistGenes[iGene]
            handAx = plt.subplot(arrayGridSpecTopHist[iGene])

            handAx.hist(np.log2(dfGeorgeSCLC[strGene].values.astype(np.float) + 1),
                        bins=np.linspace(dictGenePlotLim[strGene][0],
                                         dictGenePlotLim[strGene][1], num=31),
                        color='0.5',
                        edgecolor=None)

            handAx.set_xlim([dictGenePlotLim[strGene][0], dictGenePlotLim[strGene][1]])
            arrayXTickLoc = plt.MaxNLocator(numScatterXTicks)
            handAx.xaxis.set_major_locator(arrayXTickLoc)
            handAx.set_xticklabels([])

            handAx.set_ylim(arrayTopHistFreqLim)
            handAx.set_yticks(arrayHistFreqTicks)
            if iGene == 0:
                handAx.set_yticklabels(np.int32(arrayHistFreqTicks), fontsize=Plot.numFontSize)
                for handTick in handAx.yaxis.get_major_ticks():
                    handTick.label.set_fontsize(Plot.numFontSize)
                handAx.set_ylabel('Frequency', fontsize=Plot.numFontSize)
            else:
                handAx.set_yticklabels([])

            handAx.spines['right'].set_visible(False)
            handAx.spines['top'].set_visible(False)

            for axis in ['bottom', 'left']:
                handAx.spines[axis].set_linewidth(0.5)
            handAx.tick_params(width=0.5)

        for iGene in range(len(listSideHistGenes)):

            arrayHistFreqTicks = np.arange(start=arraySideHistFreqLim[0], stop=arraySideHistFreqLim[1], step=10)

            strGene = listSideHistGenes[iGene]
            handAx = plt.subplot(arrayGridSpecSideHist[iGene])

            handAx.hist(np.log2(dfGeorgeSCLC[strGene].values.astype(np.float) + 1),
                        bins=np.linspace(dictGenePlotLim[strGene][0],
                                         dictGenePlotLim[strGene][1], num=31),
                        color='0.5',
                        orientation='horizontal',
                        edgecolor=None)

            handAx.set_ylim([dictGenePlotLim[strGene][0], dictGenePlotLim[strGene][1]])
            arrayYTickLoc = plt.MaxNLocator(numScatterYTicks)
            handAx.yaxis.set_major_locator(arrayYTickLoc)
            handAx.set_yticklabels([])

            handAx.set_xlim(arraySideHistFreqLim)
            handAx.set_xticks(arrayHistFreqTicks)
            if iGene == 1:
                handAx.set_xticklabels(np.int32(arrayHistFreqTicks), fontsize=Plot.numFontSize)
                for handTick in handAx.xaxis.get_major_ticks():
                    handTick.label.set_fontsize(Plot.numFontSize)
                handAx.set_xlabel('Frequency', fontsize=Plot.numFontSize)
            else:
                handAx.set_xticklabels([])

            handAx.spines['right'].set_visible(False)
            handAx.spines['top'].set_visible(False)

            for axis in ['bottom', 'left']:
                handAx.spines[axis].set_linewidth(0.5)
            handAx.tick_params(width=0.5)

        handFig.savefig(os.path.join(PathDir.strOutputLoc, 'SuppFig1_-_MarkerScatter.png'),
                        ext='png', dpi=300)
        handFig.savefig(os.path.join(PathDir.strOutputLoc, 'SuppFig1_-_MarkerScatter.pdf'),
                        ext='pdf', dpi=300)
        plt.close(handFig)

        return flagResult

#   #   #   #   #
# Control the execution of various functions as required


# _ = Load.george_2015_metadata()
# _ = Load.george_2015_transcript()
# _ = PreProcess.george_2015_small_cell()

# _ = PreProcess.george_subsets()
# _ = PreProcess.george_stage()

# dictLocalNKSig = GeneSetScoring.CuratedList.nk_cells_in_tumour()


_ = Plot.figure_one()
_ = Plot.supp_fig_one_scatter()
# _ = Plot.figure_two()