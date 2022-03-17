#include <Mesh.hpp>
static Mesh * createMesh()
{
std::vector<MVertex *> vertices;
vertices.push_back(new MVertex(0, -0.0, 0.0, 0.0));
vertices.push_back(new MVertex(1, 1.0, 0.0, 0.0));
vertices.push_back(new MVertex(2, 1.0, 1.0, 0.0));
vertices.push_back(new MVertex(3, 1.0, 1.0, 1.0));
vertices.push_back(new MVertex(4, 1.0, 0.0, 1.0));
vertices.push_back(new MVertex(5, 0.0, 1.0, 0.0));
vertices.push_back(new MVertex(6, 0.0, 1.0, 1.0));
vertices.push_back(new MVertex(7, 0.0, 0.0, 1.0));
vertices.push_back(new MVertex(8, 0.5, 0.5, 0.0));
vertices.push_back(new MVertex(9, 0.5, 0.0, 0.5));
vertices.push_back(new MVertex(10, 0.0, 0.5, 0.5));
vertices.push_back(new MVertex(11, 0.5, 1.0, 0.5));
vertices.push_back(new MVertex(12, 0.5, 0.5, 1.0));
vertices.push_back(new MVertex(13, 1.0, 0.5, 0.5));
std::vector<MElement *> elements;
std::vector<MVertex *> vertices0 = {vertices[4],vertices[3],};
elements.push_back(new MElement(0, 1, vertices0));	std::vector<MVertex *> vertices1 = {vertices[0],vertices[5],vertices[10],};
elements.push_back(new MElement(1, 2, vertices1));	std::vector<MVertex *> vertices2 = {vertices[5],vertices[6],vertices[10],};
elements.push_back(new MElement(2, 2, vertices2));	std::vector<MVertex *> vertices3 = {vertices[0],vertices[10],vertices[7],};
elements.push_back(new MElement(3, 2, vertices3));	std::vector<MVertex *> vertices4 = {vertices[6],vertices[7],vertices[10],};
elements.push_back(new MElement(4, 2, vertices4));	std::vector<MVertex *> vertices5 = {vertices[12],vertices[10],vertices[13],vertices[9],};
elements.push_back(new MElement(5, 4, vertices5));	std::vector<MVertex *> vertices6 = {vertices[12],vertices[10],vertices[11],vertices[13],};
elements.push_back(new MElement(6, 4, vertices6));	std::vector<MVertex *> vertices7 = {vertices[10],vertices[13],vertices[9],vertices[8],};
elements.push_back(new MElement(7, 4, vertices7));	std::vector<MVertex *> vertices8 = {vertices[10],vertices[11],vertices[13],vertices[8],};
elements.push_back(new MElement(8, 4, vertices8));	std::vector<MVertex *> vertices9 = {vertices[1],vertices[0],vertices[9],vertices[8],};
elements.push_back(new MElement(9, 4, vertices9));	std::vector<MVertex *> vertices10 = {vertices[6],vertices[3],vertices[12],vertices[11],};
elements.push_back(new MElement(10, 4, vertices10));	std::vector<MVertex *> vertices11 = {vertices[5],vertices[2],vertices[11],vertices[8],};
elements.push_back(new MElement(11, 4, vertices11));	std::vector<MVertex *> vertices12 = {vertices[1],vertices[4],vertices[13],vertices[9],};
elements.push_back(new MElement(12, 4, vertices12));	std::vector<MVertex *> vertices13 = {vertices[7],vertices[0],vertices[10],vertices[9],};
elements.push_back(new MElement(13, 4, vertices13));	std::vector<MVertex *> vertices14 = {vertices[0],vertices[5],vertices[10],vertices[8],};
elements.push_back(new MElement(14, 4, vertices14));	std::vector<MVertex *> vertices15 = {vertices[6],vertices[7],vertices[10],vertices[12],};
elements.push_back(new MElement(15, 4, vertices15));	std::vector<MVertex *> vertices16 = {vertices[3],vertices[4],vertices[12],vertices[13],};
elements.push_back(new MElement(16, 4, vertices16));	std::vector<MVertex *> vertices17 = {vertices[2],vertices[1],vertices[13],vertices[8],};
elements.push_back(new MElement(17, 4, vertices17));	std::vector<MVertex *> vertices18 = {vertices[5],vertices[6],vertices[10],vertices[11],};
elements.push_back(new MElement(18, 4, vertices18));	std::vector<MVertex *> vertices19 = {vertices[4],vertices[7],vertices[12],vertices[9],};
elements.push_back(new MElement(19, 4, vertices19));	std::vector<MVertex *> vertices20 = {vertices[2],vertices[3],vertices[11],vertices[13],};
elements.push_back(new MElement(20, 4, vertices20));	std::vector<MVertex *> vertices21 = {vertices[13],vertices[1],vertices[9],vertices[8],};
elements.push_back(new MElement(21, 4, vertices21));	std::vector<MVertex *> vertices22 = {vertices[4],vertices[12],vertices[13],vertices[9],};
elements.push_back(new MElement(22, 4, vertices22));	std::vector<MVertex *> vertices23 = {vertices[11],vertices[2],vertices[13],vertices[8],};
elements.push_back(new MElement(23, 4, vertices23));	std::vector<MVertex *> vertices24 = {vertices[10],vertices[6],vertices[12],vertices[11],};
elements.push_back(new MElement(24, 4, vertices24));	std::vector<MVertex *> vertices25 = {vertices[7],vertices[10],vertices[12],vertices[9],};
elements.push_back(new MElement(25, 4, vertices25));	std::vector<MVertex *> vertices26 = {vertices[0],vertices[10],vertices[9],vertices[8],};
elements.push_back(new MElement(26, 4, vertices26));	std::vector<MVertex *> vertices27 = {vertices[10],vertices[5],vertices[11],vertices[8],};
elements.push_back(new MElement(27, 4, vertices27));	std::vector<MVertex *> vertices28 = {vertices[3],vertices[12],vertices[11],vertices[13],};
elements.push_back(new MElement(28, 4, vertices28));	Mesh * m =    new Mesh(3, vertices, elements);
return m;
 }; 
 