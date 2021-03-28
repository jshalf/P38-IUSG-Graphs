#ifndef MY_GLOBALS_H
#define MY_GLOBALS_H

#include<vector>

extern double totalExecutionTime;
extern double convert2SpNodeTime;
extern double computePageRankTime;
extern double prioritizeMoveTime;
extern double prioritizeSpNodeMoveTime;
extern double updateMembersTime;
extern double networkReadTime;
extern double initiationTime;
extern double calibrationTime;
extern double prioritizeMoveNodeAccessTime;
extern double prioritizeMoveNodeAccessTimeChunk;
extern double prioritizeSpNodeModuleAccessTimeChunk;
extern double bestSpModuleSelectionTime;
extern double prioritizeMoveModuleAccessTime;
extern double prioritizeSpMoveNodeAccessTime;
extern double prioritizeSpMoveModuleAccessTime;
extern double randomNumberGeneratorTime;
extern double bestModuleSelectionTime;
extern std::vector<double> threadtimes;
extern double prioritizeSpNodeCommunicationTime;
extern double prioritizeNodeCommunicationTime;
extern double prioritizeSpNodeSyncTime;
extern double prioritizeNodeSyncTime;
extern double convertModulePart1;
extern double convertModulePart2;
extern double convertModulePart3;
extern double calibrationTime;

extern int iterations_updateMembers;
extern int iteration_convertModules;
extern int iteration_prioritizeMove;
extern int iteration_prioritizeSpMove;


#endif