#ifndef ZMESH_MESH_ATTR_H
#define ZMESH_MESH_ATTR_H
#include <string>
#include <vector>
#include "Color.h"

enum AttrRate { AR_UNIFORM, AR_FACE, AR_EDGE, AR_VERTEX, AR_FACE_VERT };

enum AttrType {
    AT_DBL, AT_INT, AT_STRING, AT_VEC3, 
    AT_VEC_DBL, AT_VEC_INT, AT_VEC_VEC3, AT_VEC_BOOL, AT_VEC_COLOR, AT_VEC_LINE,
    AT_FEATURES, AT_UNKNOWN
};

class MeshAttrBase
{
public:
    MeshAttrBase(const std::string& attrName, AttrRate at_rate, AttrType at_type = AT_UNKNOWN) : mAttrName(attrName), mRate(at_rate), mType(at_type) {}
	virtual ~MeshAttrBase() = 0 {};
	virtual MeshAttrBase* clone() const = 0;

    AttrRate attrRate() const { return mRate;  }
	AttrType attrType() const { return mType; }
	const std::string& attrName() const { return mAttrName; }
	std::string& attrName() { return mAttrName; }

protected:
	std::string mAttrName;
	AttrType mType;
    AttrRate mRate;
};

template<typename T>
class MeshAttr : public MeshAttrBase
{
public:
	MeshAttr(const std::string& attrName, AttrRate rate, AttrType type = UNKNOWN_TYPE) 
        : MeshAttrBase(attrName, rate, type) {}
	MeshAttr(const T& data, const std::string& attrName, AttrRate rate, AttrType at_type = UNKNOWN_TYPE)
        : MeshAttrBase(attrName, rate, at_type), mData(data) {}
	virtual ~MeshAttr(){}
	virtual MeshAttrBase* clone() const {
		MeshAttrBase* pattr = new MeshAttr<T>(mData, mAttrName, mRate, mType);
		return pattr;
	}

	const T& attrValue() const { return mData; }
	T& attrValue() { return mData; }

protected:
	T mData;
};

class MeshFeature
{
public:
	int m_index;
	int m_scale;
	int m_note;
	double m_scalar1;
	MeshFeature() : m_index(-1), m_scale(-1), m_note(0) {}
	MeshFeature(int i) : m_index(i), m_scale(0), m_note(0) {}
	MeshFeature(int i, int s) : m_index(i), m_scale(s), m_note(0) {}
	virtual ~MeshFeature() {}
};

class MeshFeatureList
{
public:
	MeshFeatureList() {}
	MeshFeatureList(const MeshFeatureList& ml2)
	{
		this->m_vFeatures.clear();
		for (MeshFeature* pf : ml2.m_vFeatures)
			this->m_vFeatures.push_back(new MeshFeature(*pf));
	}
	~MeshFeatureList() { clear(); }
	void clear()
	{
		for (MeshFeature* f : m_vFeatures) delete f;
		m_vFeatures.clear();
	}
	int size() const { return (int)m_vFeatures.size(); }
	void addFeature(MeshFeature* mf) { m_vFeatures.push_back(mf); }
	void addFeature(int index, int scale) { m_vFeatures.push_back(new MeshFeature(index, scale)); }
	MeshFeature* back() { return m_vFeatures.back(); }
	std::vector<MeshFeature*>& getFeatureVector() { return m_vFeatures; }
	const std::vector<MeshFeature*>& getFeatureVector() const { return m_vFeatures; }
	
public:
	std::vector<MeshFeature*> m_vFeatures;
};

class LineSegment
{
public:
    LineSegment() : directional(false), color1(ZGeom::ColorBlue), color2(ZGeom::ColorRed) {}
    LineSegment(const ZGeom::Vec3d& p1, const ZGeom::Vec3d& p2, bool dir = false)
        : first(p1), second(p2), directional(dir), color1(ZGeom::ColorBlue), color2(ZGeom::ColorRed) {}
    void setColor(const ZGeom::Colorf& c) { color1 = c; }

public:    
    ZGeom::Vec3d first, second;
    bool directional;
    ZGeom::Colorf color1, color2;
};

typedef std::vector<LineSegment> MeshLineList;

typedef MeshAttr<std::vector<ZGeom::Colorf>> AttrVertColors;
typedef MeshAttr<std::vector<double>> AttrVertScalars;
typedef MeshAttr<double> AttrMeshScalar;
typedef MeshAttr<int> AttrMeshInt;
typedef MeshAttr<MeshFeatureList> AttrMeshFeatures;
typedef MeshAttr<MeshLineList> AttrMeshLines;

#endif