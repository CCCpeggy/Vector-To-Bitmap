#include "../Include/Common.h"
#pragma once
class MyObj {
private:
	cv::Mat image;
	std::vector<int> offsets;
	std::vector<CVSystem::MyPoint> vertexList;				// The actual vertices
	std::vector<CVSystem::MyIndexedTriangle> indexedTriangles;		// All triangles
	std::vector<CVSystem::MyIndexedTriangle> borderSCTriangles;		// Additional screentone triangles
	std::vector<CVSystem::MyIndexedTriangle> borderWTriangles;		// Additional white triangles
	std::vector<CVSystem::MyIndexedTriangle> scTriangles;			// Screentone triangles
	std::vector<CVSystem::MyIndexedTriangle> wTriangles;				// White triangles
	std::vector<CVSystem::MyIndexedTriangle> bTriangles;				// Black triangles
	std::vector<CVSystem::MyIndexedBezierCurves> indexedBezierCurves;	// Bezier curves (start and end points only)
	std::vector<CVSystem::MyQuad> lsQuads;
	GLuint vao;
	GLuint vbo;
	GLuint blackEbo;
	GLuint whiteEbo;
	GLuint whiteEbo2;
	GLuint scEbo;
	GLuint scEbo2;
public:
	int img_width;
	int img_height;
	MyObj(): offsets(9) {
	}

	void reset() {

	}

	void LoadFile(std::string myObjFilePath) {

		/*string myObjFileName = myObjFilePath.substr(myObjFilePath.find_last_of("/") + 1);
		string fileName = myObjFileName.substr(0, myObjFileName.find("_B.myobj"));
		string mFilePath = myObjFilePath.substr(0, myObjFilePath.find_last_of("/"));
		string filePath = mFilePath.substr(0, myObjFilePath.find_last_of("/"));
		string bFilePath = filePath  + "/Mesh/" + fileName + "_B.png";
		string cFrfilePath = filePath + "/Cartoon Filter Region/" + fileName + "_CFR.png";*/

		//check file
		if (!Common::FileExists(myObjFilePath.c_str()))
		{
			std::cout << "File: " << myObjFilePath << " is not Exist." << std::endl;
			return;
		}
		std::cout << "Loading " << myObjFilePath << "......" << std::endl;

		//if (!Common::FileExists(bFilePath.c_str()))
		//{
		//	cout << "File: " << bFilePath << " is not Exist." << endl;
		//	return;
		//}

		//if (!Common::FileExists(cFrfilePath.c_str()))
		//{
		//	cout << "File: " << cFrfilePath << " is not Exist." << endl;
		//	return;
		//}

		//int last_ = mainFilename.find_last_of("_") - 1;

		//string orifilePath;
		//string oriFile = mFile;
		//oriFile.replace(mFile.find("_CFR.", last_), 5, ".");
		//int _dpi = oriFile.find_last_of("_");
		//int _sub = oriFile.find_last_of(".");
		//oriFile.replace(_dpi, _sub - _dpi + 1, ".");
		//orifilePath = SystemParams::str_Resources_Original + "/" + oriFile;
		//if (!FileExists(orifilePath.c_str()))
		//{
		//	cout << "File: " << oriFile << " is not Exist.FF" << endl;
		//	return;
		//}
		
		ReadFromMYOBJ(myObjFilePath);
		std::cout << "Load " << myObjFilePath << " Finish." << std::endl;
		InitBuffer();
	}

	static void ReadFromMYOBJ(std::string meshFilename,
		int& img_width,
		int& img_height,
		cv::Mat& image,
		std::vector<int>& offsets,
		std::vector<CVSystem::MyPoint>& vertexList,				// The actual vertices
		std::vector<CVSystem::MyIndexedTriangle>& indexedTriangles,		// All triangles
		std::vector<CVSystem::MyIndexedTriangle>& borderSCTriangles,		// Additional screentone triangles
		std::vector<CVSystem::MyIndexedTriangle>& borderWTriangles,		// Additional white triangles
		std::vector<CVSystem::MyIndexedTriangle>& scTriangles,			// Screentone triangles
		std::vector<CVSystem::MyIndexedTriangle>& wTriangles,				// White triangles
		std::vector<CVSystem::MyIndexedTriangle>& bTriangles,				// Black triangles
		std::vector<CVSystem::MyIndexedBezierCurves>& indexedBezierCurves,	// Bezier curves (start and end points only)
		std::vector<CVSystem::MyQuad>& lsQuads)
	{
		vertexList.clear();
		indexedTriangles.clear();
		borderSCTriangles.clear();
		borderWTriangles.clear();
		scTriangles.clear();
		wTriangles.clear();
		bTriangles.clear();
		indexedBezierCurves.clear();
		lsQuads.clear();
		int nowRow = 0;

		//write binary file
		FILE* ost = fopen((meshFilename + "b").c_str(), "wb");
		bool OPB = false;
		std::ifstream myfile(meshFilename);
		while (!myfile.eof())
		{
			std::string line;
			std::getline(myfile, line);

			if (line.size() == 0) { continue; }

			std::vector<std::string> arrayStr = Common::split(line, ' ');

			// vertex list
			if (StartWith("v2d", line) && arrayStr.size() == 3)
			{
				float x = std::stof(arrayStr[1]);
				float y = std::stof(arrayStr[2]);
				vertexList.push_back(CVSystem::MyPoint(x, y));

				//write binary file
				if (OPB)
				{
					std::string str = "v2d ";

					fwrite(&str, str.length(), 1, ost);
					fwrite(&x, sizeof(float), 1, ost);
					fwrite(&y, sizeof(float), 1, ost);
				}
			}

			// all triangles
			else if (StartWith("ft", line) && arrayStr.size() == 5)
			{
				CVSystem::MyIndexedTriangle tr(std::stoi(arrayStr[1]),
					std::stoi(arrayStr[2]),
					std::stoi(arrayStr[3]),
					(CVSystem::TriangleType)std::stoi(arrayStr[4]));

				indexedTriangles.push_back(tr);

				if (tr.tri_type == CVSystem::TRIANGLE_SCREENTONE) { scTriangles.push_back(tr); }
				else if (tr.tri_type == CVSystem::TRIANGLE_WHITE) { wTriangles.push_back(tr); }
				else if (tr.tri_type == CVSystem::TRIANGLE_BLACK) { bTriangles.push_back(tr); }

				//write binary file
				if (OPB)
				{
					std::string str = "ft ";

					fwrite(&str, str.length(), 1, ost);
					fwrite(&tr.idx0, sizeof(int), 1, ost);
					fwrite(&tr.idx1, sizeof(int), 1, ost);
					fwrite(&tr.idx2, sizeof(int), 1, ost);
					fwrite(&tr.tri_type, sizeof(int), 1, ost);
				}
			}

			// border screentone triangle
			else if (StartWith("fbsc", line) && arrayStr.size() == 4)
			{
				CVSystem::MyIndexedTriangle tr(std::stoi(arrayStr[1]),
					std::stoi(arrayStr[2]),
					std::stoi(arrayStr[3]),
					CVSystem::TRIANGLE_SCREENTONE);
				borderSCTriangles.push_back(tr);

				//write binary file
				if (OPB)
				{
					std::string str = "fbsc ";

					fwrite(&str, str.length(), 1, ost);
					fwrite(&tr.idx0, sizeof(int), 1, ost);
					fwrite(&tr.idx1, sizeof(int), 1, ost);
					fwrite(&tr.idx2, sizeof(int), 1, ost);
				}
			}

			else if (StartWith("fbw", line) && arrayStr.size() == 4)
			{
				CVSystem::MyIndexedTriangle tr(std::stoi(arrayStr[1]),
					std::stoi(arrayStr[2]),
					std::stoi(arrayStr[3]),
					CVSystem::TRIANGLE_WHITE);
				borderWTriangles.push_back(tr);

				//write binary file
				if (OPB)
				{
					std::string str = "fbw ";

					fwrite(&str, str.length(), 1, ost);
					fwrite(&tr.idx0, sizeof(int), 1, ost);
					fwrite(&tr.idx1, sizeof(int), 1, ost);
					fwrite(&tr.idx2, sizeof(int), 1, ost);
				}
			}

			else if (StartWith("ibz", line))
			{
				CVSystem::MyIndexedBezierCurves bzc;
				for (size_t a = 1; a < arrayStr.size(); a++) { bzc.indices.push_back(std::stoi(arrayStr[a])); }

				//bzc.p0_points.resize(bzc.indices.size());
				//bzc.p1_points.resize(bzc.indices.size());
				//bzc.p2_points.resize(bzc.indices.size());

				indexedBezierCurves.push_back(bzc);

				//write binary file
				if (OPB)
				{
					std::string str = "ibz ";

					fwrite(&str, str.length(), 1, ost);
					for (size_t a = 1; a < arrayStr.size(); a++) { int ii = std::stoi(arrayStr[a]);  fwrite(&ii, sizeof(int), 1, ost); }
				}
			}

			else if (StartWith("q", line) && arrayStr.size() == 9)
			{
				CVSystem::MyQuad qd(CVSystem::MyPoint(std::stod(arrayStr[1]), std::stod(arrayStr[2])),
					CVSystem::MyPoint(std::stod(arrayStr[3]), std::stod(arrayStr[4])),
					CVSystem::MyPoint(std::stod(arrayStr[5]), std::stod(arrayStr[6])),
					CVSystem::MyPoint(std::stod(arrayStr[7]), std::stod(arrayStr[8])));

				lsQuads.push_back(qd);
				//lsQuads.push_back(qd);

				//write binary file
				if (OPB)
				{
					std::string str = "q ";

					fwrite(&str, str.length(), 1, ost);
					for (size_t a = 1; a < arrayStr.size(); a++) { float dd = std::stof(arrayStr[a]);  fwrite(&dd, sizeof(float), 1, ost); }
				}
			}

			else if (StartWith("iw", line) && arrayStr.size() == 2)
			{
				img_width = std::stoi(arrayStr[1]);

				//write binary file
				if (OPB)
				{
					std::string str = "iw ";

					fwrite(&str, str.length(), 1, ost);
					unsigned short int iw = img_width;
					fwrite(&iw, sizeof(iw), 1, ost);
				}
			}

			else if (StartWith("ih", line) && arrayStr.size() == 2)
			{
				img_height = std::stoi(arrayStr[1]);

				//write binary file
				if (OPB)
				{
					std::string str = "ih ";

					fwrite(&str, str.length(), 1, ost);
					unsigned short int ih = img_height;
					fwrite(&ih, sizeof(ih), 1, ost);
				}
			}

			else if (StartWith("of", line) && arrayStr.size() == 9)
			{
				offsets[0] = std::stoi(arrayStr[1]);
				offsets[1] = std::stoi(arrayStr[2]);
				offsets[2] = std::stoi(arrayStr[3]);
				offsets[3] = std::stoi(arrayStr[4]);
				offsets[4] = std::stoi(arrayStr[5]);
				offsets[5] = std::stoi(arrayStr[6]);
				offsets[6] = std::stoi(arrayStr[7]);
				offsets[7] = std::stoi(arrayStr[8]);
				//part_offset = std::stoi(arrayStr[1]);

				//write binary file
				if (OPB)
				{
					std::string str = "of ";

					fwrite(&str, str.length(), 1, ost);
					for (size_t a = 1; a < arrayStr.size(); a++) { int ii = std::stoi(arrayStr[a]);  fwrite(&ii, sizeof(ii), 1, ost); }
				}
			}

			else if (StartWith("img", line))
			{


				if (nowRow == 0)
					image = cv::Mat::zeros(img_height, img_width, CV_8UC1);
				int parserOffect = 1;
				for (int w = 0; w < img_width;)
				{
					int nowC; //0>B, 1>W
					if (arrayStr[parserOffect][0] == 'B')
					{
						nowC = 0;
					}
					else if (arrayStr[parserOffect][0] == 'W')
					{
						nowC = 255;
					}
					unsigned short int count = std::stoi(arrayStr[parserOffect + 1]);

					for (int pc = 0; pc < count; pc++)
						image.at<uchar>(nowRow, w + pc) = nowC;
					w += count;
					parserOffect += 2;
				}
				nowRow++;

			}
		}
	}

	void InitBuffer() {
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);

		std::vector<float> vertics = GetVertics(vertexList);
		glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertics.size(), &vertics[0], GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), 0);

		//glBufferData(GL_ARRAY_BUFFER, sizeof(double) * vertexList.size(), &vertexList[0], GL_STATIC_DRAW);

		//glEnableVertexAttribArray(0);
		//glVertexAttribPointer(0, 2, GL_DOUBLE, GL_FALSE, sizeof(CVSystem::MyPoint), (void*)((int)&(vertexList[0].x) - (int)&vertexList[0]));

		std::vector<unsigned int> triIndices = GetIndices(bTriangles);
		glGenBuffers(1, &blackEbo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, blackEbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * triIndices.size(), &triIndices[0], GL_STATIC_DRAW);

		triIndices = GetIndices(wTriangles);
		glGenBuffers(1, &whiteEbo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, whiteEbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * triIndices.size(), &triIndices[0], GL_STATIC_DRAW);

		triIndices = GetIndices(scTriangles);
		glGenBuffers(1, &scEbo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, scEbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * triIndices.size(), &triIndices[0], GL_STATIC_DRAW);

		triIndices = GetIndices(borderWTriangles);
		glGenBuffers(1, &whiteEbo2);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, whiteEbo2);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * triIndices.size(), &triIndices[0], GL_STATIC_DRAW);

		triIndices = GetIndices(borderSCTriangles);
		glGenBuffers(1, &scEbo2);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, scEbo2);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * triIndices.size(), &triIndices[0], GL_STATIC_DRAW);
	}

	void DrawBlack() {
		glBindVertexArray(vao);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, blackEbo);
		glDrawElements(GL_TRIANGLES, bTriangles.size() * 3, GL_UNSIGNED_INT, 0);
	}

	void DrawWhite() {
		glBindVertexArray(vao);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, whiteEbo);
		glDrawElements(GL_TRIANGLES, wTriangles.size() * 3, GL_UNSIGNED_INT, 0);
	}

	void DrawWhite2() {
		glBindVertexArray(vao);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, whiteEbo2);
		glDrawElements(GL_TRIANGLES, wTriangles.size() * 3, GL_UNSIGNED_INT, 0);
	}

	void DrawSC() {
		glBindVertexArray(vao);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, scEbo);
		glDrawElements(GL_TRIANGLES, scTriangles.size() * 3, GL_UNSIGNED_INT, 0);
	}

	void DrawSC2() {
		glBindVertexArray(vao);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, scEbo2);
		glDrawElements(GL_TRIANGLES, scTriangles.size() * 3, GL_UNSIGNED_INT, 0);
	}

private:

	static std::vector<unsigned int> GetIndices(std::vector<CVSystem::MyIndexedTriangle>& myIndexedTriangle) {
		std::vector<unsigned int> indices;
		for (auto tri : myIndexedTriangle) {
			indices.push_back(tri.idx0);
			indices.push_back(tri.idx1);
			indices.push_back(tri.idx2);
		}
		return indices;
	}

	static std::vector<float> GetVertics(std::vector<CVSystem::MyPoint>& pointList) {
		std::vector<float> vertics;
		for (auto p : pointList) {
			vertics.push_back(p.x);
			vertics.push_back(p.y);
		}
		return vertics;
	}

	void ReadFromMYOBJ(std::string meshFilename)
	{
		ReadFromMYOBJ(meshFilename, img_width, img_height, image, offsets, vertexList, indexedTriangles,
			borderSCTriangles, borderWTriangles, scTriangles, wTriangles, bTriangles, indexedBezierCurves, lsQuads);
	}

	// string startwith
	static bool StartWith(std::string prefix, std::string argument)
	{
		if (argument.substr(0, prefix.size()) == prefix)
		{
			return true;
		}
		return false;
	}
};