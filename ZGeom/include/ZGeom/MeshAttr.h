#ifndef ZMESH_MESH_ATTR_H
#define ZMESH_MESH_ATTR_H
#include <string>
#include <iostream>

enum AttrRate {UNIFORM, FACE, EDGE, VERTEX, FACE_VERT};
enum AttrType {CPP_DOUBLE, CPP_INT, CPP_STRING, CPP_VECTOR_DOUBLE, CPP_VECTOR_INT, CPP_COLOR, UNKNOWN_TYPE};

class MeshAttrBase
{
public:
	MeshAttrBase(AttrRate rate, const std::string& attrName, AttrType type = UNKNOWN_TYPE) : mRate(rate), mAttrName(attrName), mType(type) {}
	virtual ~MeshAttrBase() = 0 {};
	const std::string& getAttrName() const { return mAttrName; }
	void setAttrName(const std::string& attrName) { mAttrName = attrName; }
	AttrRate attrRate() const { return mRate; }
	AttrType attrType() const { return mType; }
	virtual MeshAttrBase* clone() const = 0;

protected:
	std::string mAttrName;
	AttrRate mRate;
	AttrType mType;
};

template<typename T>
class MeshAttr : public MeshAttrBase
{
public:
	MeshAttr(AttrRate rate, const std::string& attrName, AttrType type = UNKNOWN_TYPE) :  MeshAttrBase(rate, attrName, type){}
	MeshAttr(const T& data, AttrRate rate, const std::string& attrName, AttrType type = UNKNOWN_TYPE) : MeshAttrBase(rate, attrName, type), mData(data) {}
	virtual ~MeshAttr(){}
	virtual MeshAttrBase* clone() const {
		MeshAttrBase* pattr = new MeshAttr<T>(mData, mRate, mAttrName, mType);
		return pattr;
	}

	T getValue() const { return mData; }
	T& getValue() { return mData; }

private:
	T mData;
};

#endif