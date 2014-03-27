#ifndef ZMESH_MESH_ATTR_H
#define ZMESH_MESH_ATTR_H
#include <string>
#include <vector>
#include "Color.h"

enum AttrRate {AR_UNIFORM, AR_FACE, AR_EDGE, AR_VERTEX, AR_FACE_VERT};
enum AttrType {AT_DBL, AT_INT, AT_STRING, AT_VEC_DBL, AT_VEC_INT, AT_VEC_COLOR, AT_FEATURES, AT_UNKNOWN};

class MeshAttrBase
{
public:
	MeshAttrBase(const std::string& attrName, AttrType at_type = AT_UNKNOWN) : mAttrName(attrName), mType(at_type) {}
	virtual ~MeshAttrBase() = 0 {};
	virtual MeshAttrBase* clone() const = 0;

	AttrType attrType() const { return mType; }
	const std::string& attrName() const { return mAttrName; }
	std::string& attrName() { return mAttrName; }

protected:
	std::string mAttrName;
	AttrType mType;
};

template<typename T, AttrRate R>
class MeshAttr : public MeshAttrBase
{
public:
	MeshAttr(const std::string& attrName, AttrType type = UNKNOWN_TYPE) :  MeshAttrBase(attrName, type){}
	MeshAttr(const T& data, const std::string& attrName, AttrType at_type = UNKNOWN_TYPE) : MeshAttrBase(attrName, at_type), mData(data) {}
	virtual ~MeshAttr(){}
	virtual MeshAttrBase* clone() const {
		MeshAttrBase* pattr = new MeshAttr<T,R>(mData, mAttrName, mType);
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
	MeshFeature() : m_index(-1), m_scale(-1), m_note(0) {}
	MeshFeature(int i) : m_index(i), m_scale(0), m_note(0) {}
	MeshFeature(int i, int s) : m_index(i), m_scale(s), m_note(0) {}
	virtual ~MeshFeature() {}
};

class MeshFeatureList
{
public:
	~MeshFeatureList() { clear(); }

	void addFeature(MeshFeature* mf) { m_vFeatures.push_back(mf); }
	std::vector<MeshFeature*>* getFeatureVector() { return &m_vFeatures; }

	void clear() 
	{
		for (MeshFeature* f : m_vFeatures) delete f;
		m_vFeatures.clear();
	}

	int featureType;
	std::vector<MeshFeature*> m_vFeatures;
};

typedef MeshAttr<std::vector<ZGeom::Colorf>, AR_VERTEX> AttrVertColors;
typedef MeshAttr<std::vector<double>, AR_VERTEX> AttrVertScalars;
typedef MeshAttr<double, AR_UNIFORM> AttrMeshScalar;
typedef MeshAttr<int, AR_UNIFORM> AttrMeshInt;
typedef MeshAttr<MeshFeatureList, AR_UNIFORM> AttrMeshFeatures;

#endif