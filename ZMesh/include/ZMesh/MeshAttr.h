#ifndef ZMESH_MESH_ATTR_H
#define ZMESH_MESH_ATTR_H
#include <string>
#include <iostream>

enum AttrRate {UNIFORM, FACE, EDGE, VERTEX};
class MeshAttrBase
{
public:
    MeshAttrBase(AttrRate rate, const std::string& attrName) : mRate(rate), mAttrName(attrName) {}
    virtual ~MeshAttrBase() = 0 { std::cout << "Destroy MeshAttrBase!" << std::endl; };
    const std::string& getAttrName() const { return mAttrName; }
    void setAttrName(const std::string& attrName) { mAttrName = attrName; }

protected:
    std::string mAttrName;
    AttrRate mRate;
};

template<typename T>
class MeshAttr : public MeshAttrBase
{
public:
    MeshAttr(AttrRate rate, const std::string& attrName) :  MeshAttrBase(rate, attrName){}
    MeshAttr(const T& data, AttrRate rate, const std::string& attrName) : MeshAttrBase(rate, attrName), mData(data) {}
    virtual ~MeshAttr(){ std::cout << "Destroy MeshAttr!" << std::endl; }

    T getValue() const { return mData; }
    T& getValue() { return mData; }

private:
    T mData;
};

#endif