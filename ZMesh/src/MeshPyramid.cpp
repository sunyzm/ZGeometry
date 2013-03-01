#include <cmath>
#include <ctime>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <exception>
#include <fstream>
#include "MeshPyramid.h"

#define SILVER_TRIANGLE_TRHESH 0.05
#define AVOID_SILVER_TRIANGLE
#define DEBUG_OUTPUT
#define BOUNDARY_PENALTY	1 ///1000
#define COST_MIN			1e15	

using namespace std;
	

MeshPyramid::MeshPyramid( CMesh* mesh )
{
	setInitialMesh(mesh);
}

MeshPyramid::MeshPyramid()
{
	originalMesh = NULL;
	m_levels = 0;
	m_Id2IndexMap = NULL;
}

void MeshPyramid::setInitialMesh( CMesh* mesh )
{
	if(!mesh->m_bIsPointerVectorExist)
		throw logic_error("Pointer Vector not exist for MeshPyramid::setInitialMesh()");

	originalMesh = mesh;

	m_nVertices = originalMesh->m_nVertex;
	m_nFace = originalMesh->m_nFace;

	//calculate initial quadric for each face and vertex
	m_FaceQuadrics.clear();
	m_VertexQuadrics.clear();
	m_FaceQuadrics.resize(m_nFace);
	m_VertexQuadrics.reserve(m_nVertices);

	for (int i = 0; i < m_nFace; ++i)
	{
		CFace* pFace = originalMesh->m_vFaces[i];
		vector<double> funcPara = pFace->getPlaneFunction();
		Quadric q(funcPara[0], funcPara[1], funcPara[2], funcPara[3], pFace->getArea());
		q *= q.getArea();
		m_FaceQuadrics[i] = q;
	}

	// calculate initial quadric for each vertex
	for (int i = 0; i < m_nVertices; ++i)
	{
		const CVertex* pV = originalMesh->m_vVertices[i];
		// iterate through each face incident on v; 
		// TODO: consider boundary constraint
		vector<CFace*> adjFaces = pV->getAdjacentFaces();
		Quadric q;
		for (unsigned int j = 0; j < adjFaces.size(); ++j)
			q += m_FaceQuadrics[adjFaces[j]->m_fIndex];
		m_VertexQuadrics.push_back(q);

	}

	//add boundary constraints on every boundary edge
	for (int i = 0; i < originalMesh->m_nHalfEdge; ++i)
	{
		const CHalfEdge* pHE = originalMesh->m_vHalfEdges[i];
		if (pHE->m_eTwin && pHE->m_eTwin->m_bIsValid)
			continue;

		int bV1 = pHE->m_Vertices[0]->getVID(), bv2 = pHE->m_Vertices[1]->getVID();

		Vector3D vf = pHE->m_Face->getNormal();
		Vector3D ve = pHE->m_Vertices[1]->getPosition() - pHE->m_Vertices[0]->getPosition();
		Vector3D vn = vf ^ ve;
		vn.normalize();
		double para_a = vn[0], 
			   para_b = vn[1],
			   para_c = vn[2],
			   para_d = vn * pHE->m_Vertices[0]->getPosition(),
			   para_area = ve.length2(); 
		Quadric bQ = Quadric(para_a, para_b, para_c, para_d, para_area);
		bQ *= bQ.getArea() * BOUNDARY_PENALTY;
				
// 		int faceIndex = pHE->m_Face->m_fIndex;
// 		Quadric bQ = m_FaceQuadrics[faceIndex];

		m_VertexQuadrics[bV1-1] += bQ;
		m_VertexQuadrics[bv2-1] += bQ;
	}
}

MeshPyramid::~MeshPyramid()
{
	clear();
}

void MeshPyramid::buildPyramid()
{
	if (!originalMesh) throw logic_error("Initial mesh required for building pyramid!");

	/* construct multi-resolution mesh structure */
	{
		MeshLevel baseLevel;
		baseLevel.mesh = originalMesh;
		m_vMeshes.clear();
		m_vMeshes.push_back(baseLevel);
	}
	
	for (int level = 1; level < m_levels; ++level)
	{
		MeshLevel meshLevel;
		meshLevel.mesh = new CMesh;
		CMesh* currentMesh = meshLevel.mesh;
		currentMesh->clone(m_vMeshes[level-1].mesh);
		
		list<VertexPair> vertexPairsList;
		//initialize vertex pairs
		{
			bool *isHalfEdgeVisited = new bool[currentMesh->m_nEdge];
			for (int i = 0; i < currentMesh->m_nEdge; ++i)
				isHalfEdgeVisited[i] = false;
			for (int i = 0; i < currentMesh->m_nEdge; ++i)
			{
				CHalfEdge* curHE = currentMesh->m_HalfEdges[i];
				if (!isHalfEdgeVisited[i])
				{
					VertexPair vp;
					vp.pV[0] = currentMesh->m_HalfEdges[i]->m_Vertices[0];
					vp.pV[1] = currentMesh->m_HalfEdges[i]->m_Vertices[1];
					vp.pHE = currentMesh->m_HalfEdges[i];
					Quadric q = m_VertexQuadrics[vp.pV[0]->m_vid-1] + m_VertexQuadrics[vp.pV[1]->m_vid-1];
					double cost_0 = q.eval(vp.pV[0]->m_vPosition),	//cost of pV[1] collapse to pV[0]
						   cost_1 = q.eval(vp.pV[1]->m_vPosition);	//cost of pV[0] collapse to pV[1]
					if (cost_0 < cost_1)
					{
						vp.vToKeep = 0;
						vp.cost = cost_0;
					} else
					{
						vp.vToKeep = 1;
						vp.cost = cost_1;
					}
					vertexPairsList.push_back(vp);

					if (currentMesh->m_HalfEdges[i]->m_eTwin && currentMesh->m_HalfEdges[i]->m_eTwin->m_bIsValid)
						isHalfEdgeVisited[currentMesh->m_HalfEdges[i]->m_eTwin->m_eIndex] = true;
				}
			}
			delete []isHalfEdgeVisited;
		}

		if (level == 1) ////build initial verticesCoverFull		
		{			
			for (int i = 1; i <= currentMesh->m_nVertex; ++i)
			{
				//assign initial vertex cover; each vertex only cover itself
				list<int> coverList;
				coverList.push_back(i);		//here i is fixed vid, not variable vIndex
				meshLevel.verticesCoverFull.insert(make_pair(i, coverList));
			}
		}	//if level == 1
		else meshLevel.verticesCoverFull = m_vMeshes[level-1].verticesCoverFull;

		/************* simplify current mesh by half *****************/
		/*iteration:
			select least cost pair
			do pair contraction
			update cost for remaining pair */
		
		/* find initial min-cost pair */
		double minCost = COST_MIN;		//need more consideration
		list<VertexPair>::iterator minIter;
		for (list<VertexPair>::iterator vp_iter = vertexPairsList.begin(); vp_iter != vertexPairsList.end(); ++vp_iter)
		{
			if (vp_iter->cost < minCost)
			{
				minCost = vp_iter->cost;
				minIter = vp_iter;
			}
		}
		
		/** main iteration for each edge collapse **/
//		const int num_contraction = int(this->m_nVertices / pow(2.0, level)); 
		const int num_contraction = int(3 * this->m_nVertices / pow(4.0, level)); 
		for (int step = 0; step < num_contraction; ++step)		
		{
			CVertex* vKeep = minIter->pV[minIter->vToKeep];
			CVertex* vRemove = minIter->pV[1 - minIter->vToKeep];
			const VID idKeep = vKeep->m_vid; 
			const VID idRemove = vRemove->m_vid;
			m_VertexQuadrics[idKeep-1] += m_VertexQuadrics[idRemove-1];		//combine quadric of vKeep and vRemove
			
			CHalfEdge *heRemove1(NULL), *heRemove2(NULL);
			CFace *fRemove1(NULL), *fRemove2(NULL);
			CHalfEdge* eDelIn[3] = {NULL, NULL, NULL}; 
			CHalfEdge* eDelOut[3] = {NULL, NULL, NULL};
			CFace *fIn(NULL), *fOut(NULL);
			CVertex *vOppoIn(NULL), *vOppoOut(NULL);

			heRemove1 = minIter->correspondingHalfEdge();
			assert(heRemove1);
			fRemove1 = heRemove1->m_Face;
			if (heRemove1->m_eTwin && heRemove1->m_eTwin->m_bIsValid)
			{
				heRemove2 = heRemove1->m_eTwin;
				fRemove2 = heRemove2->m_Face;
			}
			
			bool bVertexPairValid(true);
			if (/*!heRemove1 || */ !minIter->pV[0]->m_bIsValid || !minIter->pV[1]->m_bIsValid)
			{
				cout << "Invalid vertex pair!" << endl;
				bVertexPairValid = false;
			}
			if (!currentMesh->isHalfEdgeMergeable(heRemove1))
			{
				cout << "Half-edge not mergeable!" << endl;
				bVertexPairValid = false;
			}

#ifdef AVOID_SILVER_TRIANGLE
			/* judge whether silver triangle will appear */
			else {
				for (vector<CHalfEdge*>::iterator he_iter = vRemove->m_HalfEdges.begin(); he_iter != vRemove->m_HalfEdges.end(); ++he_iter)
				{
					CFace* curF = (*he_iter)->m_Face;
					if (curF == fRemove1 || curF == fRemove2)
						continue;
					CVertex *v1 = (*he_iter)->m_eNext->m_Vertices[0],
							*v2 = (*he_iter)->m_eNext->m_Vertices[1],
							*v3 = vKeep;
					Vector3D l1 = v1->getPos() - v2->getPos(),
							 l2 = v2->getPos() - v3->getPos(),
							 l3 = v3->getPos() - v1->getPos();
					double area = (l1 ^ l2).length() / 2;
					double gamma = 4.0 * std::sqrt(3.0) * area / (l1.length2() + l2.length2() + l3.length2());
					if (gamma < SILVER_TRIANGLE_TRHESH)
					{
						//cout << "Silver Triangle " << gamma << " avoided!" << endl;
						bVertexPairValid = false;
						break;
					}
				}
			}
#endif	//avoid silver triangle
			if (!bVertexPairValid)
			{
				vertexPairsList.erase(minIter);
				minCost = COST_MIN;
				for (list<VertexPair>::iterator vp_iter = vertexPairsList.begin(); 
					vp_iter != vertexPairsList.end(); ++vp_iter)
				{
					if (vp_iter->cost < minCost)
					{
						minCost = vp_iter->cost;
						minIter = vp_iter;
					}
				}
				continue;
			}
			
			const bool bInRemove = (heRemove1->m_Vertices[1] == vRemove);
			if (bInRemove)
			{
				eDelIn[0] = heRemove1; 
				eDelIn[1] = heRemove1->m_eNext; 
				eDelIn[2] = heRemove1->m_ePrev;
				fIn = eDelIn[0]->m_Face;
				vOppoIn = eDelIn[1]->m_Vertices[1];
				if (heRemove2)
				{
					eDelOut[0] = heRemove2; 
					eDelOut[1] = heRemove2->m_eNext; 
					eDelOut[2] = heRemove2->m_ePrev;
					fOut = eDelOut[0]->m_Face;
					vOppoOut = eDelOut[1]->m_Vertices[1];
				}
			}
			else 
			{
				eDelOut[0] = heRemove1; 
				eDelOut[1] = heRemove1->m_eNext; 
				eDelOut[2] = heRemove1->m_ePrev;
				fOut = eDelOut[0]->m_Face;
				vOppoOut = eDelOut[1]->m_Vertices[1];
				if (heRemove2)
				{
					eDelIn[0] = heRemove2; 
					eDelIn[1] = heRemove2->m_eNext; 
					eDelIn[2] = heRemove2->m_ePrev;
					fIn = eDelIn[0]->m_Face;
					vOppoIn = eDelIn[1]->m_Vertices[1];
				}
			}
						
			assert(vRemove->m_HalfEdges.size() > 0);

			//case 1
			if (vRemove->m_HalfEdges.size() == 1)
			{
				CHalfEdge* eDel[3] = {vRemove->m_HalfEdges[0], vRemove->m_HalfEdges[0]->m_eNext, vRemove->m_HalfEdges[0]->m_ePrev};
				CVertex *v1 = vRemove->m_HalfEdges[0]->m_Vertices[1],
					    *v2 = vRemove->m_HalfEdges[0]->m_ePrev->m_Vertices[0];
				/* remove obsolete half-edges from vertex */
				{
					vector<CHalfEdge*>::iterator he_iter1 = find(v1->m_HalfEdges.begin(), v1->m_HalfEdges.end(), eDel[1]);
					vector<CHalfEdge*>::iterator he_iter2 = find(v2->m_HalfEdges.begin(), v2->m_HalfEdges.end(), eDel[2]);
					v1->m_HalfEdges.erase(he_iter1);
					v2->m_HalfEdges.erase(he_iter2);
				}

				if (eDel[1]->m_eTwin && eDel[1]->m_eTwin->m_bIsValid)
					eDel[1]->m_eTwin->m_eTwin = NULL;

				fRemove1->m_bIsValid = false;
				for (int k = 0; k < 3; ++k)
				{
					eDel[k]->m_bIsValid = false;
				}
				
			} 
			//case 2
			else /* if (vRemove->m_HalfEdges.size() >= 2)*/
			{
				if (eDelIn[0] && eDelIn[0]->m_bIsValid)
				{
					vector<CHalfEdge*>::iterator he_iter1 = find(vOppoIn->m_HalfEdges.begin(), vOppoIn->m_HalfEdges.end(), eDelIn[2]);
					eDelIn[1]->m_Vertices[1]->m_HalfEdges.erase(he_iter1);
					vector<CHalfEdge*>::iterator he_iter3 = find(vKeep->m_HalfEdges.begin(), vKeep->m_HalfEdges.end(), eDelIn[0]);
					vKeep->m_HalfEdges.erase(he_iter3);
				}
				if (eDelOut[0] && eDelOut[0]->m_bIsValid)
				{
					vector<CHalfEdge*>::iterator he_iter2 = find(vOppoOut->m_HalfEdges.begin(), vOppoOut->m_HalfEdges.end(), eDelOut[2]);
					eDelOut[1]->m_Vertices[1]->m_HalfEdges.erase(he_iter2);
					vector<CHalfEdge*>::iterator he_iter4 = find(vKeep->m_HalfEdges.begin(), vKeep->m_HalfEdges.end(), eDelOut[1]);
					vKeep->m_HalfEdges.erase(he_iter4);
				}

				for (vector<CHalfEdge*>::iterator iter = vRemove->m_HalfEdges.begin(); iter != vRemove->m_HalfEdges.end(); ++iter)
				{
					CHalfEdge* curHe = *iter;
					assert(curHe && curHe->m_bIsValid);
					CFace* curF = curHe->m_Face;
					if (curF == fIn)
					{
						//assert(curHe == eDelIn[1]);
						continue;
					}
					if (curF == fOut)
					{
						//assert(curHe == eDelOut[0]);
						continue;
					}
					
					curHe->m_Vertices[0] = curHe->m_ePrev->m_Vertices[1] = vKeep;
					replace(curF->m_Vertices.begin(), curF->m_Vertices.end(), vRemove, vKeep);

					vKeep->m_HalfEdges.push_back(curHe);
				}

				if (fIn)
				{
					if (eDelIn[2]->m_eTwin/* && eDelIn[2]->m_eTwin->m_bIsValid*/)
						eDelIn[2]->m_eTwin->m_eTwin = eDelIn[1]->m_eTwin;

					if (eDelIn[1]->m_eTwin/* && eDelIn[1]->m_eTwin->m_bIsValid*/)
						eDelIn[1]->m_eTwin->m_eTwin = eDelIn[2]->m_eTwin;
				}				
				if (fOut)
				{
					if (eDelOut[1]->m_eTwin/* && eDelOut[1]->m_eTwin->m_bIsValid*/)
						eDelOut[1]->m_eTwin->m_eTwin = eDelOut[2]->m_eTwin;

					if (eDelOut[2]->m_eTwin/* && eDelOut[2]->m_eTwin->m_bIsValid*/)
						eDelOut[2]->m_eTwin->m_eTwin = eDelOut[1]->m_eTwin;
				}

				//invalidate edges and faces
				if (fIn)
				{
					eDelIn[0]->m_bIsValid = eDelIn[1]->m_bIsValid = eDelIn[2]->m_bIsValid = false;
					fIn->m_bIsValid = false;
				}
				if (fOut)
				{
					eDelOut[0]->m_bIsValid = eDelOut[1]->m_bIsValid = eDelOut[2]->m_bIsValid = false;
					fOut->m_bIsValid = false;
				}
			}//end of case 2
				
			/*** update verticesCover ***/			
			VerticesCover::iterator iterKeepFull = meshLevel.verticesCoverFull.find(idKeep),
 					                iterRemoveFull = meshLevel.verticesCoverFull.find(idRemove);
 			iterKeepFull->second.insert(iterKeepFull->second.end(), iterRemoveFull->second.begin(), iterRemoveFull->second.end());
 			meshLevel.verticesCoverFull.erase(iterRemoveFull);
			
			VerticesCover::iterator iterKeep = meshLevel.verticesCover.find(idKeep);
			if (iterKeep == meshLevel.verticesCover.end())
			{
				list<int> clist;
				clist.push_back(idRemove);
				meshLevel.verticesCover.insert(make_pair(idKeep, clist));
				iterKeep = meshLevel.verticesCover.find(idKeep);
			}
			else iterKeep->second.push_back(idRemove);
			
			vRemove->m_bIsValid = false;

			if (vOppoIn && vOppoIn->m_HalfEdges.size() == 0)
			{
				iterRemoveFull = meshLevel.verticesCoverFull.find(vOppoIn->m_vid);
				iterKeepFull->second.insert(iterKeepFull->second.end(), iterRemoveFull->second.begin(), iterRemoveFull->second.end());
				meshLevel.verticesCoverFull.erase(iterRemoveFull);
				
				iterKeep->second.push_back(vOppoIn->m_vid);
				vOppoIn->m_bIsValid = false;
			}
			if (vOppoOut && vOppoOut->m_HalfEdges.size() == 0)
			{
				iterRemoveFull = meshLevel.verticesCoverFull.find(vOppoOut->m_vid);
				iterKeepFull->second.insert(iterKeepFull->second.end(), iterRemoveFull->second.begin(), iterRemoveFull->second.end());
				meshLevel.verticesCoverFull.erase(iterRemoveFull);
				
				iterKeep->second.push_back(vOppoOut->m_vid);
				vOppoOut->m_bIsValid = false;
			}

			//** update vertexPairsList **//
			vertexPairsList.erase(minIter);
			minCost = COST_MIN;
			for (list<VertexPair>::iterator vp_iter = vertexPairsList.begin(); 
				vp_iter != vertexPairsList.end(); )
			{
				if (vp_iter->pV[0] == vRemove)
				{
					if (vp_iter->pV[1] == vOppoIn || vp_iter->pV[1] == vOppoOut || vp_iter->pV[1] == vKeep)
					{
						vp_iter = vertexPairsList.erase(vp_iter);
						continue;
					}
					else
					{
						vp_iter->pV[0] = vKeep;
						vp_iter->pHE = vp_iter->correspondingHalfEdge();
						assert(vp_iter->pHE);
					}
				}
				if (vp_iter->pV[1] == vRemove)
				{
					if (vp_iter->pV[0] == vOppoIn || vp_iter->pV[0] == vOppoOut || vp_iter->pV[0] == vKeep)
					{
						vp_iter = vertexPairsList.erase(vp_iter);
						continue;
					}
					else
					{
						vp_iter->pV[1] = vKeep;
						vp_iter->pHE = vp_iter->correspondingHalfEdge();
						assert(vp_iter->pHE);
					}
				}
				//assert(vp_iter->pV[0] != vp_iter->pV[1]);
				if ( !vp_iter->pV[0]->m_bIsValid || !vp_iter->pV[1]->m_bIsValid )
				{
					cout << "Vertex pair has void vertex!" << endl;
					vp_iter = vertexPairsList.erase(vp_iter);
					continue;
				}
								
				if (vp_iter->pV[0] == vKeep || vp_iter->pV[1] == vKeep)
				{
					Quadric q = m_VertexQuadrics[vp_iter->pV[0]->m_vid - 1] + m_VertexQuadrics[vp_iter->pV[1]->m_vid - 1];
					double cost_0 = q.eval(vp_iter->pV[0]->m_vPosition),
						cost_1 = q.eval(vp_iter->pV[1]->m_vPosition);
					if (cost_0 < cost_1)
					{
						vp_iter->vToKeep = 0;
						vp_iter->cost = cost_0;
					} else
					{
						vp_iter->vToKeep = 1;
						vp_iter->cost = cost_1;
					}
				}
				if (vp_iter->cost < minCost)
				{
					minCost = vp_iter->cost;
					minIter = vp_iter;
				}
				vp_iter++;
			}

		}////for every contraction

		/////////// after simplification by half ////////////////////////////////////////////////
		////strip out invalid vertices, half-edges and faces; re-assign index
		for (vector<CHalfEdge*>::iterator iter = currentMesh->m_HalfEdges.begin(); iter != currentMesh->m_HalfEdges.end();iter++)
		{
			CHalfEdge* curHe = *iter;
			if (curHe->m_eTwin && !curHe->m_eTwin->m_bIsValid)
				curHe->m_eTwin = NULL;
		}
		for (vector<CHalfEdge*>::iterator iter = currentMesh->m_HalfEdges.begin(); iter != currentMesh->m_HalfEdges.end();)
		{
			CHalfEdge* curHe = *iter;
			if (!curHe->m_bIsValid)
			{
				delete curHe;
				iter = currentMesh->m_HalfEdges.erase(iter);
				continue;
			}
			assert(curHe->m_Vertices[0]->m_bIsValid && curHe->m_Vertices[1]->m_bIsValid);
			iter++;
		}
		currentMesh->m_nEdge = (int)currentMesh->m_HalfEdges.size();
		
		for (vector<CFace*>::iterator iter = currentMesh->m_Faces.begin(); iter != currentMesh->m_Faces.end();)
		{
			CFace* curF = *iter;
			if (curF->m_bIsValid == false)
			{
				delete curF;
				iter = currentMesh->m_Faces.erase(iter);
				continue;
			}
			if(curF->m_Vertices.empty())
			{
				assert(curF->m_HalfEdges.empty());
				cout << "Empty face!" << endl;
				delete curF;
				iter = currentMesh->m_Faces.erase(iter);
				continue;
			}
			for (int j = 0; j < 3; ++j)
			{
				assert(curF->m_Vertices[j]->m_bIsValid);
			}
			iter++;
		}
		currentMesh->m_nFace = (int)currentMesh->m_Faces.size();

		for (vector<CVertex*>::iterator iter = currentMesh->m_Vertices.begin(); iter != currentMesh->m_Vertices.end();)
		{
			CVertex* curV = *iter;
			if (curV->m_bIsValid == false)
			{
				delete curV;
				iter = currentMesh->m_Vertices.erase(iter);
			}
			else if (curV->m_HalfEdges.size() == 0)
			{
				cout << "Isolated vertex: " << curV->m_vid << endl;
				delete curV;
				iter = currentMesh->m_Vertices.erase(iter);
			}
			else
			{
				curV->m_nValence = int(curV->m_HalfEdges.size());
				curV->judgeOnBoundary();
				iter++;
			}
		}
		currentMesh->m_nVertex = (int)currentMesh->m_Vertices.size();
		currentMesh->m_nFace = (int)currentMesh->m_Faces.size();
		
		cout << "  Vertex Num: " << currentMesh->m_nVertex << endl
			//<< "  Face Num: " << currentMesh->m_nFace << endl
			//<< "  Half-edge Num: " << currentMesh->m_nEdge << endl
			//<< "  Edge Num: " << currentMesh->GetEdgeNum() << endl 
			<< "  Boundary Vertex Num: " << currentMesh->getBoundaryVertexNum() << endl; 
		
		currentMesh->assignElementsIndex();
		currentMesh->buildArrayRepresentation();
		currentMesh->gatherStatistics();

		//finally
		m_vMeshes.push_back(meshLevel);
	}///for each level; outer most iteration

	/*** build id2index map ***/
	m_Id2IndexMap = new int*[m_nVertices];
	for (int i = 0; i < m_nVertices; ++i)
	{
		m_Id2IndexMap[i] = new int[m_levels];
		for (int l = 1; l < m_levels; ++l)
			m_Id2IndexMap[i][l] = -1;
	}
	
	for (int i = 0; i < m_nVertices; ++i)
		m_Id2IndexMap[i][0] = i;

	for (int l = 1; l < m_levels; ++l)
	{
		const CMesh* mesh_l = m_vMeshes[l].mesh;
		const VerticesCover& vc = m_vMeshes[l].verticesCoverFull;
		for (map<VID, list<VID> >::const_iterator miter = vc.begin(); miter != vc.end(); ++miter)
		{
			int cvid = miter->first;
			const list<VID>& vlist = miter->second;

			int matchedIdx = -1;
			for (int j = 0; j < mesh_l->m_nVertex; ++j)
			{
			 	if (mesh_l->m_pVertex[j].m_vid == cvid)
			 	{
			 		matchedIdx = j;
			 		break;
			 	}
			}
			if (matchedIdx == -1)
			{
				cout << "No Matched Index in level " << l << "!!!" << endl;
				exit(-1);
			}

			for (list<VID>::const_iterator liter = vlist.begin(); liter != vlist.end(); ++liter)
			{
				int vid = *liter;
				m_Id2IndexMap[vid-1][l] = matchedIdx;
			}
		}
	}
}

void MeshPyramid::clear()
{
	for (int i = 0; i < m_nVertices; ++i)
		delete []m_Id2IndexMap[i];
	delete []m_Id2IndexMap;

	for (unsigned int i = 1; i < m_vMeshes.size(); ++i)	//leave original mesh for outer destruction
	{
		if (m_vMeshes[i].mesh != NULL) 
			delete m_vMeshes[i].mesh;
	}
	for (unsigned int i = 0; i < m_vProlongMat.size()-1; ++i)
	{
		if (m_vProlongMat[i]) delete m_vProlongMat[i];
	}
	m_vMeshes.clear();
	m_FaceQuadrics.clear();
	m_VertexQuadrics.clear();
	m_vProlongMat.clear();
	originalMesh = NULL;
	m_levels = 0;
}

CMesh* MeshPyramid::getMesh( int level ) const
{
	assert(m_levels == (int)m_vMeshes.size());
	assert( -m_levels <= level && level < m_levels);
	if (level >= 0) 
		return m_vMeshes[level].mesh;
	else 
		return m_vMeshes[level + m_levels].mesh;
}

void MeshPyramid::buildProlongMatrix()
{
	clock_t start, end;
	start = clock();

	for (int n = 0; n < m_vMeshes.size()-1; ++n)
	{
		SparseMatrix *smat = new SparseMatrix();
		// the mat should be v_n * v_(n+1)

		const CMesh* mesh_n = m_vMeshes[n].mesh, *mesh_n1 = m_vMeshes[n+1].mesh;  
		for (int i = 0; i < mesh_n->m_nVertex; ++i)
		{
			const CVertex* vi = mesh_n->m_Vertices[i];
			double minDistToFace = 1e15;
			CFace*  minFace(NULL);
			vector<double> minDistCoord;
			vector<double> tmpCoord;
			tmpCoord.resize(3);

			//shortcut to vertex that has direct correspondence; seems working
			int vid = vi->getVID();
			bool hasCorrespondence(false);
			for (vector<CVertex*>::const_iterator viter = mesh_n1->m_Vertices.begin(); viter != mesh_n1->m_Vertices.end(); ++viter)
			{
				if (vid == (*viter)->getVID())
				{
					smat->m_ii.push_back(i + 1);
					smat->m_jj.push_back((*viter)->m_vIndex + 1);
					smat->m_ss.push_back(1.0);
					hasCorrespondence = true;
					break;
				}
			}
			if (hasCorrespondence) continue;

			for (int j = 0; j < mesh_n1->m_nFace; ++j)
			{
				CFace* fj = mesh_n1->m_Faces[j];
				double dist = fj->distanceToVertex(vi, tmpCoord);
				if (dist < minDistToFace)
				{
					minFace = fj;
					minDistToFace = dist;
					minDistCoord = tmpCoord;
				}
			}
			int faceIndex[3] = {
						minFace->m_Vertices[0]->m_vIndex, 
						minFace->m_Vertices[1]->m_vIndex,
						minFace->m_Vertices[2]->m_vIndex
						};
			for (int l = 0; l < 3; ++l)
			{
				if ( minDistCoord[l] < 1e-6) continue;	//zero element ignored to save space
				smat->m_ii.push_back(i+1);
				smat->m_jj.push_back(faceIndex[l] + 1);
				smat->m_ss.push_back(minDistCoord[l]);
			}
			
		}//	for (int i = 0; i < mesh_n->m_nVertex; ++i)
		m_vProlongMat.push_back(smat);
	}

	end = clock();
	double duration = double(end-start) / CLOCKS_PER_SEC;
	cout << "Prolongation matrix build time: " << duration << endl;
}

SparseMatrix* MeshPyramid::getProlongationMatrix( int level ) const
{
	assert (0 <= level && level <= m_levels - 2 && m_vProlongMat.size() == m_levels - 1);
	return m_vProlongMat[level];
}

CHalfEdge* MeshPyramid::VertexPair::correspondingHalfEdge()
{
	CHalfEdge* ret(NULL);
	for(vector<CHalfEdge*>::iterator iter = pV[0]->m_HalfEdges.begin(); 
		iter != pV[0]->m_HalfEdges.end();
		++iter)
	{
		if ((*iter)->m_Vertices[1] == pV[1] && (*iter)->m_bIsValid)
		{
			ret = *iter;
			return ret;
		}
	}
	if (ret == NULL)
	{
		for (vector<CHalfEdge*>::iterator iter = pV[1]->m_HalfEdges.begin();
			 iter != pV[1]->m_HalfEdges.end(); ++iter)
		{
			if ( (*iter)->m_Vertices[1] == pV[0] && (*iter)->m_bIsValid )
			{
				ret = *iter;
				return ret;
			}
		}
	}
	
	return NULL;
}

void SparseMatrix::dump( std::string filename ) const
{
	assert(m_ii.size() == m_jj.size() && m_ii.size() == m_ss.size());
	ofstream spmOut(filename.c_str(), ios::trunc);
	for (size_t i = 0; i < m_ii.size(); ++i)
		spmOut << m_ii[i] << '\t' << m_jj[i] << '\t' << m_ss[i] << endl;
	spmOut.close();
}


void MeshPyramid::dumpVertexValence( int level, std::string filename )
{
	ofstream valOut(filename.c_str(), ios::trunc);
	CMesh* tmesh = getMesh(level);

	for (int i = 0; i < tmesh->m_nVertex; ++i)
	{
		valOut << i << '\t' << tmesh->m_Vertices[i]->m_nValence << '\t' << tmesh->m_pVertex[i].m_nValence << endl;
	}

	valOut.close();
}

std::list<int> MeshPyramid::getCoveredVertexList( int level, int idx ) const
{
	assert( 0 <= level && level < m_levels-1);

	const VerticesCover& coverMap = m_vMeshes[level+1].verticesCover;
	const CMesh* tmesh = getMesh(level);
	VID v_id = tmesh->m_pVertex[idx].m_vid;
	map<VID, list<VID> >::const_iterator viter = coverMap.find(v_id);
	
	list<VID> retlist;
	if (viter == coverMap.end()) return retlist;

	const list<VID>& vlist = viter->second;
	for (list<VID>::const_iterator iter = vlist.begin(); iter != vlist.end(); ++iter)
	{
		int fromID = *iter;
		int toIdx = m_Id2IndexMap[fromID-1][level];
		retlist.push_back(toIdx);
	}

	return retlist;
}
