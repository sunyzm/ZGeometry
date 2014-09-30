#include "Geodesic.h"
#include <vector>
using namespace std;

class GeoNote
{
public:
    int m_id;
    double m_geodesic;
public:
    GeoNote(int mid, double geo) { m_id = mid; m_geodesic = geo; }
    GeoNote& operator = (const GeoNote& note) { m_id = note.m_id; m_geodesic = note.m_geodesic; return(*this); }
    friend bool operator > (const GeoNote& note1, const GeoNote& note2) { return note1.m_geodesic > note2.m_geodesic; }
};
typedef std::priority_queue<GeoNote, std::vector<GeoNote>, std::greater<GeoNote> > GeoQueue;

struct VertHelper{
    int m_mark;
    double m_LocalGeodesic;
    bool m_inheap;
    
    VertHelper() : m_mark(-1), m_LocalGeodesic(-1.), m_inheap(false) {}
};

double calLocalGeodesic(const CMesh& mesh, const vector<VertHelper>& vVertNotes, int ia, int ib, int ic)
{   
    // ia - vertex with smaller geodesic; ib - with greater geodesic; ic - update
    double la = (mesh.getVertexPosition(ib) - mesh.getVertexPosition(ic)).length();
    double lb = (mesh.getVertexPosition(ia) - mesh.getVertexPosition(ic)).length();
    double lc = (mesh.getVertexPosition(ib) - mesh.getVertexPosition(ia)).length();
    double ctheta = (la*la + lb*lb - lc*lc) / (2 * la*lb);
    double stheta = sqrt(1 - ctheta*ctheta);
    double u = vVertNotes[ib].m_LocalGeodesic - vVertNotes[ia].m_LocalGeodesic;
    double ld = lb - la*ctheta;
    double le = la*stheta;
    double delta = lc*lc - u*u;
    double tc = vVertNotes[ic].m_LocalGeodesic;
    if (delta >= 0.0) {
        delta = sqrt(delta);
        double t1 = lb*(u*ld + le*delta) / (lc*lc);//(-B+delta)/(2*A);
        if (t1 > u && lb*(t1 - u) / t1 > la*ctheta && lb*(t1 - u) / t1 < la / abs(ctheta))
        {
            if (tc < 0) tc = t1 + vVertNotes[ia].m_LocalGeodesic;
            else tc = min(tc, t1 + vVertNotes[ia].m_LocalGeodesic);
        }
        else {
            double minab = min(lb + vVertNotes[ia].m_LocalGeodesic, la + vVertNotes[ib].m_LocalGeodesic);
            if (tc < 0) tc = minab;
            else tc = min(tc, minab);
        }
    }
    else {
        double minab = min(lb + vVertNotes[ia].m_LocalGeodesic, la + vVertNotes[ib].m_LocalGeodesic);
        if (tc < 0) tc = minab;
        else tc = min(tc, minab);
    }

    return tc;
}

double ZGeom::calGeodesic(const CMesh& mesh, int s, int t)
{
    assert(s >= 0 && t >= 0 && s < mesh.vertCount() && t < mesh.vertCount());
    if (s == t) return 0.0;
    const CVertex& notei = *mesh.getVertex(s);
    GeoQueue heapqueue; 
    vector<VertHelper> vVertNotes(mesh.vertCount());    
    vVertNotes[s].m_mark = s;
    vVertNotes[s].m_LocalGeodesic = 0;

    int j, k, ia, ib, ic;
    double geo = 0.0;

    for (j = 0; j < notei.outValence(); ++j) {
        auto ee = notei.getHalfEdge(j);
        int endv = ee->getVertIndex(1);
        Vector3D vt = mesh.getVertexPosition(endv) - mesh.getVertexPosition(s);
        double mgeo = vt.length();
        if (endv == t) return mgeo; // destination reached
        vVertNotes[endv].m_LocalGeodesic = mgeo;
        vVertNotes[endv].m_mark = s;
    }

    // first ring
    for (j = 0; j < notei.outValence(); ++j) {
        auto e1 = notei.getHalfEdge(j);
        ia = e1->getVertIndex(1);
        auto e2 = e1->nextHalfEdge();
        ib = e2->getVertIndex(1);
        if (vVertNotes[ia].m_LocalGeodesic > vVertNotes[ib].m_LocalGeodesic) 
            std::swap(ia, ib);
        e1 = e2->twinHalfEdge();
        if (e1 == NULL) continue;        
        e2 = e1->nextHalfEdge();
        ic = e2->getVertIndex(1);
        double mgeo = calLocalGeodesic(mesh, vVertNotes, ia, ib, ic);
        vVertNotes[ic].m_LocalGeodesic = mgeo;
        heapqueue.push(GeoNote(ic, mgeo));
    } 

    bool stop = false;
    while (!stop && !heapqueue.empty())
    {
        GeoNote nt = heapqueue.top();
        heapqueue.pop();
        int sg = nt.m_id;
        if (sg == t) return nt.m_geodesic;

        double sgd = nt.m_geodesic;
        if (vVertNotes[sg].m_mark == s) continue; //matched already
        vVertNotes[sg].m_mark = s;

        // update adjacent vertices of sg
        const CVertex *vsg = mesh.getVertex(sg);
        for (k = 0; k < vsg->outValence(); ++k) {
            ia = sg;
            auto e1 = vsg->getHalfEdge(k);
            ib = e1->getVertIndex(1);
            if (vVertNotes[ib].m_mark != s) continue; // unreached point
            auto e2 = e1->nextHalfEdge();
            if (vVertNotes[ia].m_LocalGeodesic > vVertNotes[ib].m_LocalGeodesic) 
                std::swap(ia, ib);            
            ic = e2->getVertIndex(1);
            if (vVertNotes[ic].m_mark != s) {
                double gg = calLocalGeodesic(mesh, vVertNotes, ia, ib, ic);   // update geodesic
                if (vVertNotes[ic].m_LocalGeodesic < 0 || gg < vVertNotes[ic].m_LocalGeodesic) {
                    // heaped, shorter patch came                
                    vVertNotes[ic].m_LocalGeodesic = gg;
                    heapqueue.push(GeoNote(ic, gg));
                }
            }
            e2 = e1->twinHalfEdge();
            if (!e2) continue;
            e1 = e2->nextHalfEdge();
            ic = e1->getVertIndex(1);
            if (vVertNotes[ic].m_mark != s) {
                double gg = calLocalGeodesic(mesh, vVertNotes, ia, ib, ic); // update geodesic
                if (vVertNotes[ic].m_LocalGeodesic < 0 || gg < vVertNotes[ic].m_LocalGeodesic) {
                    // heaped, shorter patch came                
                    vVertNotes[ic].m_LocalGeodesic = gg;
                    heapqueue.push(GeoNote(ic, gg));
                }
            }
        }
    }

    return geo;
}

double ZGeom::calGeodesicToBoundary(CMesh& mesh, int s)
{
    assert(s >= 0 && s < mesh.vertCount());
    int vertCount = mesh.vertCount();
    const CVertex& notei = *mesh.getVertex(s);
    GeoQueue heapqueue;
    vector<VertHelper> vVertNotes(mesh.vertCount());
    vVertNotes[s].m_mark = s;
    vVertNotes[s].m_LocalGeodesic = 0;

    int j, k, ia, ib, ic;
    double geo = 0.0;
    auto vVertIsHole = mesh.getVertsOnHole();
    auto vVertOnBoundary = mesh.getVertsOnBoundary();
    vector<bool> vNonHoleBoundaryVert(vertCount, false);
    int nonHoleBoundary(0);
    for (int i = 0; i < vertCount; ++i) {
        if (vVertOnBoundary[i] && !vVertIsHole[i]) {
            vNonHoleBoundaryVert[i] = true;
            nonHoleBoundary++;
        }
    }    
    if (0 == nonHoleBoundary) return 0;

    for (j = 0; j < notei.outValence(); ++j) {
        auto ee = notei.getHalfEdge(j);
        int endv = ee->getVertIndex(1);
        Vector3D vt = mesh.getVertexPosition(endv) - mesh.getVertexPosition(s);
        double mgeo = vt.length();
        if (vNonHoleBoundaryVert[endv]) return mgeo;
        vVertNotes[endv].m_LocalGeodesic = mgeo;
        vVertNotes[endv].m_mark = s;
    }

    for (j = 0; j < notei.outValence(); ++j) {
        auto e1 = notei.getHalfEdge(j);
        ia = e1->getVertIndex(1);
        auto e2 = e1->nextHalfEdge();
        ib = e2->getVertIndex(1);
        if (vVertNotes[ia].m_LocalGeodesic > vVertNotes[ib].m_LocalGeodesic)
            swap(ia, ib);
        e1 = e2->twinHalfEdge();
        if (e1 == NULL) continue;
        e2 = e1->nextHalfEdge();
        ic = e2->getVertIndex(1);
        double mgeo = calLocalGeodesic(mesh, vVertNotes, ia, ib, ic);
        vVertNotes[ic].m_LocalGeodesic = mgeo;
        heapqueue.push(GeoNote(ic, mgeo));
    }

    bool stop(false);
    int count = -vertCount;
    while (!stop && !heapqueue.empty())
    {
        if (++count == 0) break;
        GeoNote nt = heapqueue.top();
        heapqueue.pop();
        int sg = nt.m_id;
        if (vNonHoleBoundaryVert[sg]) return nt.m_geodesic;

        double sgd = nt.m_geodesic;
        if (vVertNotes[sg].m_mark == s) continue;
        vVertNotes[sg].m_mark = s;
        if (vVertIsHole[sg]) continue;
        CVertex* vsg = mesh.getVertex(sg);
        for (k = 0; k < vsg->outValence(); ++k) {
            ia = sg;
            const CHalfEdge* e1 = vsg->getHalfEdge(k);
            ib = e1->getVertIndex(1);
            if (vVertNotes[ib].m_mark != s) continue;
            auto e2 = e1->nextHalfEdge();
            if (vVertNotes[ia].m_LocalGeodesic > vVertNotes[ib].m_LocalGeodesic)
                std::swap(ia, ib);
            ic = e2->getVertIndex(1);
            if (vVertNotes[ic].m_mark != s) {
                double gg = calLocalGeodesic(mesh, vVertNotes, ia, ib, ic);
                if (vVertNotes[ic].m_LocalGeodesic < 0 || gg < vVertNotes[ic].m_LocalGeodesic) {
                    vVertNotes[ic].m_LocalGeodesic = gg;
                    heapqueue.push(GeoNote(ic, gg));
                }
            }

            e2 = e1->twinHalfEdge();
            if (e2 == NULL) continue;
            e1 = e2->nextHalfEdge();
            ic = e1->getVertIndex(1);
            if (vVertNotes[ic].m_mark != s) {
                double gg = calLocalGeodesic(mesh, vVertNotes, ia, ib, ic);
                if (vVertNotes[ic].m_LocalGeodesic < 0 || gg < vVertNotes[ic].m_LocalGeodesic) {
                    vVertNotes[ic].m_LocalGeodesic = gg;
                    heapqueue.push(GeoNote(ic, gg));
                }
            }
        }        
    }

    return geo;
}

