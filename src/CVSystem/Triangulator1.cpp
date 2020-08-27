
/**
*
* Reza Adhitya Saputra (reza.adhitya.saputra@gmail.com)
* Version: 2014
*
*/

#include "../../Include/Common.h"
using namespace CVSystem;

void CVSystem::Triangulator1::TraceImage(std::string strFilename, cv::Mat img, int* mask, int* dilatedMask, bool isLog)
{
	// filename
	this->_strFilename = strFilename;

	if (isLog)
		std::cout << "Initializing Triangulation\n";

	cv::Mat invMaskImg;
	if (mask)
	{
		this->_has_tone = true;
		cv::Mat imgMaskDilated = cv::Mat::zeros(img.rows, img.cols, CV_8UC1);
		for (int a = 0; a < imgMaskDilated.cols; a++)
		{
			for (int b = 0; b < imgMaskDilated.rows; b++)
			{
				if (dilatedMask[a + b * img.cols] != -1) { imgMaskDilated.ptr<uchar>(b, a)[0] = 255; }
			}
		}

		cv::Mat imgMask = cv::Mat::zeros(img.rows, img.cols, CV_8UC1);
		for (int a = 0; a < imgMask.cols; a++)
		{
			for (int b = 0; b < imgMask.rows; b++)
			{
				if (mask[a + b * img.cols] != -1) { imgMask.ptr<uchar>(b, a)[0] = 255; }
			}
		}

		cv::bitwise_not(imgMask, invMaskImg);
		cv::distanceTransform(invMaskImg, this->_distance_img, CV_DIST_L2, CV_DIST_MASK_5);
		invMaskImg.convertTo(invMaskImg, CV_8UC1);
		this->_mask = imgMask.clone();
		this->_dilated_mask = imgMaskDilated.clone();

		img += imgMask;

		imgMask.release();
		imgMaskDilated.release();
	}
	// Dimension
	this->_w_scaled = img.cols * SystemParams::t_scale_factor;
	this->_h_scaled = img.rows * SystemParams::t_scale_factor;
	cv::Size sz(_w_scaled, _h_scaled);

	if (this->_has_tone)
	{
		this->_mask_scaled = this->_mask;
		this->_dilated_mask_scaled = this->_dilated_mask;
		this->_distance_img_scaled = this->_distance_img;
	}
	cv::Mat inpImg = img.clone();

	potrace_bitmap_t* bm = GetBM(inpImg);

	// Tracing
	potrace_param_t* param = potrace_param_default();
	param->opticurve = SystemParams::t_opticurve/*1*/;
	param->opttolerance = SystemParams::t_opttolerance/*0.2*/;
	param->progress.callback = NULL;
	_trace_state = potrace_trace(param, bm);
	//cv::imwrite("FFFF.png", inpImg);

	if (mask)
	{
		STimg = invMaskImg.clone();
		potrace_bitmap_t* bmST = GetBM(invMaskImg);
		// TracingST
		potrace_param_t* paramST = potrace_param_default();
		paramST->opticurve = SystemParams::t_opticurve;
		paramST->opttolerance = SystemParams::t_opttolerance;
		paramST->progress.callback = NULL;
		_trace_stateST = potrace_trace(paramST, bmST);

		potrace_path_t* p = _trace_state->plist;
		potrace_path_t* lastp;
		for (; p; p = p->next) { lastp = p; }
		//lastp->next = _trace_stateST->plist;
	}
	// Free resources
	potrace_param_free(param);
	free(bm->map);
	free(bm);

	inpImg.release();

	//
	//this->_has_tone = false;
	//
}

void CVSystem::Triangulator1::LSCalculate2(cv::Mat img, bool isLog, float offsetX, float offsetY)
{
	traceImg = cv::Mat(img.size(), CV_8UC3, cv::Scalar(255, 255, 255));
	traceImgS = cv::Mat(img.size(), CV_8UC3, cv::Scalar(255, 255, 255));
	// Timing calculation
	using namespace boost::chrono;
	auto total_smoothing_time = 0;
	auto ls_time = 0;

	std::vector<MyTriangle>			  scTriangles;
	std::vector<MyTriangle>			  whiteTriangles;
	std::vector<MyTriangle>			  blackTriangles;
	std::vector<std::vector<MyPoint>> allPolygons;
	std::vector<std::vector<MyPoint>> STPolygons;

	std::vector<std::vector<MyPoint>> oriPoly;
	// Get Seed
	for (potrace_path_t* p = _trace_state->plist; p; p = p->next)
	{
		if (p->sign == '-') continue;
		std::vector<MyPoint> pPolySimple;
		// Timing
		auto startGauss1 = steady_clock::now();
		std::vector<MyPoint> pPolyDense = GetPolygonLS(p);
		// Timing
		auto durGauss1 = steady_clock::now() - startGauss1;
		total_smoothing_time += duration_cast<milliseconds>(durGauss1).count();

		if (pPolyDense.size() == 0) { continue; }

		std::vector<bool> pFlag;
		RDPSimplification(pPolyDense, pPolySimple, pFlag);

		if (pPolySimple.size() > 2)
		{
			auto startLS1 = steady_clock::now();

			std::vector<MyQuad> tempQuads = GetBezierQuadLS(pPolyDense, pFlag, true);
			this->_lsQuads.insert(this->_lsQuads.begin(), tempQuads.begin(), tempQuads.end());

			auto durLS1 = steady_clock::now() - startLS1;
			ls_time += duration_cast<milliseconds>(durLS1).count();
		}
		std::vector<std::vector<MyPoint>> cPolySimple;

		for (potrace_path_t* q = p->childlist; q; q = q->sibling)
		{
			std::vector<MyPoint> qPolyTempSimple;

			// Timing
			auto start2 = steady_clock::now();

			std::vector<MyPoint> qPolyTempDense = GetPolygonLS(q);

			// Timing
			auto durGauss2 = steady_clock::now() - start2;
			total_smoothing_time += duration_cast<milliseconds>(durGauss2).count();

			if (qPolyTempDense.size() == 0) { continue; }
			std::vector<bool> qFlagTemp;
			RDPSimplification(qPolyTempDense, qPolyTempSimple, qFlagTemp);
			if (qPolyTempSimple.size() > 2)
			{
				// Timing
				auto startLS2 = steady_clock::now();

				std::vector<MyQuad> tempQuads = GetBezierQuadLS(qPolyTempDense, qFlagTemp, false);
				this->_lsQuads.insert(this->_lsQuads.begin(), tempQuads.begin(), tempQuads.end());
				if (q->sign == '+')
					cPolySimple.push_back(qPolyTempSimple);
				else
				{//reverse for ReduceSTpoint()
					std::vector<MyPoint> qPolySimpleCounterclockwise(qPolyTempSimple.rbegin(), qPolyTempSimple.rend());
					cPolySimple.push_back(qPolySimpleCounterclockwise);
				}

				// Timing
				auto durLS2 = steady_clock::now() - startLS2;
				ls_time += duration_cast<milliseconds>(durLS2).count();
			}
		}

		if (pPolySimple.size() == 0) { continue; }
		//std::cout<<"T "<<pPolySimple.size();
		//for(size_t a = 0; a < cPolySimple.size(); a++) { std::cout<<" "<<cPolySimple[a].size(); }
		CD_Cdt cdt_blk = GetCDT(pPolySimple, cPolySimple);
		//std::cout<<" CDT\n";
		allPolygons.push_back(pPolySimple);
		for (size_t a = 0; a < cPolySimple.size(); a++) { allPolygons.push_back(cPolySimple[a]); }

		for (CD_Cdt::Finite_faces_iterator fit = cdt_blk.finite_faces_begin(); fit != cdt_blk.finite_faces_end(); ++fit)
		{
			if (fit->is_in_domain())
			{
				MyTriangle tr = MyTriangle(MyPoint(fit->vertex(0)->point().x(), fit->vertex(0)->point().y()),
					MyPoint(fit->vertex(1)->point().x(), fit->vertex(1)->point().y()),
					MyPoint(fit->vertex(2)->point().x(), fit->vertex(2)->point().y()),
					TRIANGLE_NOT_DEFINED);
				int pixelValue = PixelsInTriangles::MedianValue(img, tr.p1, tr.p2, tr.p3);
				MyPoint ctr = tr.GetCenter();
				if (pixelValue < 127 || img.ptr<uchar>((int)ctr.y, (int)ctr.x)[0] < 127)
				{
					tr.tri_type = TRIANGLE_BLACK;
					blackTriangles.push_back(tr);
				}
				else {
					tr.tri_type = TRIANGLE_WHITE;
					whiteTriangles.push_back(tr);
				}
			}
		}
	}

	////////////////////
	//for debug

	cv::Mat AllStep = cv::Mat(img.size(), CV_8UC3, cv::Scalar(200, 200, 200));
	for (int a = 0; a < AllStep.cols; a++)
	{
		for (int b = 0; b < AllStep.rows; b++)
		{
			if (img.ptr<uchar>(b, a)[0] < 127)
			{
				AllStep.at<cv::Vec3b>(b, a) = cv::Vec3b(0, 0, 0);
			}
			else
			{
				AllStep.at<cv::Vec3b>(b, a) = cv::Vec3b(100, 100, 100);
			}
		}
	}

	for (int pIdx = 0; pIdx < allPolygons.size(); pIdx++)
	{
		for (int Idx = 0; Idx < allPolygons[pIdx].size() - 1; Idx++)
		{
			cv::circle(AllStep, cv::Point2d(allPolygons[pIdx][Idx].x, allPolygons[pIdx][Idx].y), 3, cv::Scalar(0, 0, 255), -1);
		}
	}
	testImg = AllStep.clone();

	////////////////////

	std::cout << "ST\n";

	if (this->_has_tone)
	{
		testcount = 0;
		for (potrace_path_t* p = _trace_stateST->plist; p; p = p->next)
		{
			if (p->sign == '-') continue;
			std::vector<MyPoint> pPolySimple;
			// Timing
			auto startGauss1 = steady_clock::now();
			std::vector<MyPoint> pPolyDense = GetPolygonLS(p);
			// Timing
			auto durGauss1 = steady_clock::now() - startGauss1;
			total_smoothing_time += duration_cast<milliseconds>(durGauss1).count();

			if (pPolyDense.size() == 0) { continue; }

			std::vector<bool> pFlag;
			RDPSimplification(pPolyDense, pPolySimple, pFlag);

			if (pPolySimple.size() > 2)
			{
				auto startLS1 = steady_clock::now();

				std::vector<MyQuad> tempQuads = GetBezierQuadLS(pPolyDense, pFlag, true);
				this->_lsQuads.insert(this->_lsQuads.begin(), tempQuads.begin(), tempQuads.end());

				auto durLS1 = steady_clock::now() - startLS1;
				ls_time += duration_cast<milliseconds>(durLS1).count();
			}
			std::vector<std::vector<MyPoint>> cPolySimple;

			for (potrace_path_t* q = p->childlist; q; q = q->sibling)
			{
				std::vector<MyPoint> qPolyTempSimple;

				// Timing
				auto start2 = steady_clock::now();

				std::vector<MyPoint> qPolyTempDense = GetPolygonLS(q);

				// Timing
				auto durGauss2 = steady_clock::now() - start2;
				total_smoothing_time += duration_cast<milliseconds>(durGauss2).count();

				if (qPolyTempDense.size() == 0) { continue; }
				std::vector<bool> qFlagTemp;
				RDPSimplification(qPolyTempDense, qPolyTempSimple, qFlagTemp);
				if (qPolyTempSimple.size() > 2)
				{
					// Timing
					auto startLS2 = steady_clock::now();

					std::vector<MyQuad> tempQuads = GetBezierQuadLS(qPolyTempDense, qFlagTemp, false);
					this->_lsQuads.insert(this->_lsQuads.begin(), tempQuads.begin(), tempQuads.end());

					if (q->sign == '+')
						cPolySimple.push_back(qPolyTempSimple);
					else
					{//reverse for ReduceSTpoint()
						std::vector<MyPoint> qPolySimpleCounterclockwise(qPolyTempSimple.rbegin(), qPolyTempSimple.rend());
						cPolySimple.push_back(qPolySimpleCounterclockwise);
					}

					// Timing
					auto durLS2 = steady_clock::now() - startLS2;
					ls_time += duration_cast<milliseconds>(durLS2).count();
				}
			}

			if (pPolySimple.size() == 0) { continue; }

			std::vector<MyPoint> rdPoly;
			ReduceSTpoint(allPolygons, pPolySimple, rdPoly);
			if (ReduceSTpoint(allPolygons, rdPoly, rdPoly) != 0)
				STPolygons.push_back(rdPoly);

			std::vector<std::vector<MyPoint>> rdcPoly;
			for (size_t a = 0; a < cPolySimple.size(); a++)
			{
				std::vector<MyPoint> stPoly;
				ReduceSTpoint(allPolygons, cPolySimple[a], stPoly);
				if (ReduceSTpoint(allPolygons, stPoly, stPoly) != 0)
					STPolygons.push_back(stPoly);
				rdcPoly.push_back(stPoly);
				stPoly.clear();
			}

			CD_Cdt cdt_blk = GetCDT(rdPoly, rdcPoly, true, true);
			rdPoly.clear();
			for (size_t a = 0; a < rdcPoly.size(); a++) { rdcPoly[a].clear(); }
			rdcPoly.clear();
			for (CD_Cdt::Finite_faces_iterator fit = cdt_blk.finite_faces_begin(); fit != cdt_blk.finite_faces_end(); ++fit)
			{
				if (fit->is_in_domain())
				{
					MyTriangle tr = MyTriangle(MyPoint(fit->vertex(0)->point().x(), fit->vertex(0)->point().y()),
						MyPoint(fit->vertex(1)->point().x(), fit->vertex(1)->point().y()),
						MyPoint(fit->vertex(2)->point().x(), fit->vertex(2)->point().y()),
						TRIANGLE_NOT_DEFINED);
					scTriangles.push_back(tr);
				}
			}
		}
	}
	//cv::imwrite("TraceImage.png", traceImg);

	//cv::imwrite("TraceImageS.png", traceImgS);
	// Timing
	if (isLog)
		std::cout << "SMoothing time: " << total_smoothing_time << " milliseconds\n";

	// Timing
	auto start = steady_clock::now();

	CD_Cdt cd_cdt1 = GetSTCDTWithRectangleBorder(allPolygons, STPolygons);	// conforming delaunay
	std::vector<CD_Cdt::Face_handle> fsh;
	if (this->_has_tone)
	{
		if (isLog)
			std::cout << "Calculate screentone borders...\n";

		// Rollback !
		for (size_t a = 0; a < fsh.size(); a++)
		{
			fsh[a]->set_marked(false);
		}
		// A seed on the border
		std::list<CD_Point> seeds;
		seeds.push_back(CD_Point(1, 1));

		CGAL::refine_Delaunay_mesh_2(cd_cdt1, seeds.begin(), seeds.end(), CD_Criteria(SystemParams::t_delaunay_aspect_bound, SystemParams::t_delaunay_max_length));

		AssignVertexIndex(cd_cdt1);
		DiscardVertexOffset(offsetX, offsetY);

		PointInTrianglesTest* bTest = new PointInTrianglesTest();
		bTest->PushTriangles(blackTriangles);

		PointInTrianglesTest* wTest = new PointInTrianglesTest();
		wTest->PushTriangles(whiteTriangles);

		PointInTrianglesTest* scTest = new PointInTrianglesTest();
		scTest->PushTriangles(scTriangles);

		for (CD_Cdt::Finite_faces_iterator fit = cd_cdt1.finite_faces_begin(); fit != cd_cdt1.finite_faces_end(); ++fit)
		{
			if (fit->is_in_domain())
			{
				TriangleType triType = TRIANGLE_NOT_DEFINED;
				if (!blackTriangles.empty() && bTest->IsInside(fit)) { triType = TRIANGLE_BLACK; }
				else if (!whiteTriangles.empty() && wTest->IsInside(fit)) { triType = TRIANGLE_WHITE; }
				else if (!scTriangles.empty() && scTest->IsInside(fit)) { triType = TRIANGLE_SCREENTONE; }
				else { triType = TRIANGLE_WHITE; }

				fit->triangleType = triType;

				MyIndexedTriangle tri(fit->vertex(0)->info(), fit->vertex(1)->info(), fit->vertex(2)->info(), triType);

				this->_indexedTriangles.push_back(tri);
				if (triType == TRIANGLE_SCREENTONE) { this->_scTriangles.push_back(tri); }
				else if (triType == TRIANGLE_WHITE) { this->_wTriangles.push_back(tri); }
				else if (triType == TRIANGLE_BLACK) { this->_bTriangles.push_back(tri); }
			}
		}
		ProcessEdgesFromCDT(cd_cdt1);
	}
	else
	{

		// A seed on the border
		std::list<CD_Point> seeds;
		seeds.push_back(CD_Point(1, 1));

		CGAL::refine_Delaunay_mesh_2(cd_cdt1, seeds.begin(), seeds.end(), CD_Criteria(SystemParams::t_delaunay_aspect_bound, SystemParams::t_delaunay_max_length));
		AssignVertexIndex(cd_cdt1);
		DiscardVertexOffset(offsetX, offsetY);

		PointInTrianglesTest* bTest = new PointInTrianglesTest();
		bTest->PushTriangles(blackTriangles);

		for (CD_Cdt::Finite_faces_iterator fit = cd_cdt1.finite_faces_begin(); fit != cd_cdt1.finite_faces_end(); ++fit)
		{
			if (fit->is_in_domain())
			{
				TriangleType triType = (bTest->IsInside(fit)) ? TRIANGLE_BLACK : TRIANGLE_WHITE;
				MyIndexedTriangle tri(fit->vertex(0)->info(), fit->vertex(1)->info(), fit->vertex(2)->info(), triType);
				this->_indexedTriangles.push_back(tri);
				fit->triangleType = triType;

				if (triType == TRIANGLE_WHITE) { this->_wTriangles.push_back(tri); }
				else if (triType == TRIANGLE_BLACK) { this->_bTriangles.push_back(tri); }
			}
		}

		ProcessEdgesFromCDT(cd_cdt1);
	}

	GetCornersFromPolygon();
	ProcessPolygonCorners();

	// Timing
	auto dur = steady_clock::now() - start;

	if (isLog)
	{
		std::cout << "least square time: " << ls_time << " milliseconds\n";
		std::cout << "triangulation time: " << duration_cast<milliseconds>(dur).count() << " milliseconds\n";
	}

	std::vector<std::string> fullpath = UtilityFunctions::split(_strFilename, '//');
	std::vector<std::string> nameArray = UtilityFunctions::split(fullpath[fullpath.size() - 1], '.');
	std::string meshFilename = SystemParams::str_Resources_Mesh + nameArray[0] + ".myobj";

	std::vector<MyQuad> resizedLsQuads = this->GetLSQuads();

	_offsets.resize(8);
	for (size_t a = 0; a < _offsets.size(); a++)
	{
		_offsets[a] = 0;
	}

	cv::Mat deOffset = cv::Mat(cv::Size(img.cols + 2 * offsetX, img.rows + 2 * offsetX), CV_8UC1, cv::Scalar(255));
	cv::Mat ROI = img(cv::Rect(std::abs(offsetX), std::abs(offsetY), deOffset.cols, deOffset.rows));
	ROI.copyTo(deOffset);
	_TriImage = deOffset.clone();
	//OBJIO::WriteToMYOBJ(meshFilename,
	//	(_w_scaled + 2 * offsetX) / SystemParams::t_scale_factor,
	//	(_h_scaled + 2 * offsetY) / SystemParams::t_scale_factor,
	//	deOffset,
	//	_offsets,
	//	_vertexList,
	//	_indexedBezierCurves,
	//	_borderSCTriangles,
	//	_borderWTriangles,
	//	resizedLsQuads,		// resized
	//	_indexedTriangles);

	//OBJWriter::WriteToOBJ("file.obj", cd_cdt1);	
	//MakeIndexedTriangles(cd_cdt1);

	if (isLog)
		std::cout << "Done\n";
}

std::vector<CVSystem::MyQuad> CVSystem::Triangulator1::GetBezierQuadLS(std::vector<CVSystem::MyPoint> poly, std::vector<bool> rdpFlags, bool isInside)
{
	std::vector<MyQuad> quads;

	// Begin to sample points
	std::vector<std::vector<MyPoint>> rdpPoly;
	rdpPoly.push_back(std::vector<CVSystem::MyPoint>());

	//MyPoint firstPt = poly[0];
	//MyPoint lastPt = poly[poly.size() - 1];

	for (size_t a = 0; a < poly.size(); a++)
	{
		int idx = rdpPoly.size() - 1;
		rdpPoly[idx].push_back(poly[a]);
		if (rdpFlags[a] && a != 0 && a != poly.size() - 1)
		{
			rdpPoly.push_back(std::vector<MyPoint>());
			rdpPoly[idx + 1].push_back(poly[a]);
		}
	}

	for (size_t a = 0; a < rdpPoly.size(); a++)
	{
		std::vector<MyPoint> segm = rdpPoly[a];

		// RESUBDIVIDE
		if (segm.size() < 5) { Resubdivide(segm); }

		std::vector<MyPoint> fourPoints;
		if (CurveFitting::PointstoBezier(segm, fourPoints))
		{
			quads.push_back(MyQuad(segm[0], fourPoints[1], fourPoints[2], segm[segm.size() - 1]));
		}
		else {}
	}

	// testing for c1 continuity
	/*MyQuad q0 = quads[0];
	MyQuad qN = quads[quads.size() - 1];

	std::cout << q0.p0.x << " " << q0.p0.y << "\n";
	std::cout << qN.p3.x << " " << qN.p3.y << "\n";
	std::cout << "\n";*/

	// Todo: enforce non-corners to be c1
	/*
	for(size_t a = 0; a < quads.size(); a++)
	{
	MyQuad* q1 = &quads[a];
	MyQuad* q2 = 0;
	if(a < quads.size() - 1) {q2 = &quads[a+1];}
	else {q2 = &quads[0];}

	// modify q1 p3
	// modify q2 p0

	MyPoint midPoint = (q1->p2 + q2->p1) / 2.0;
	q1->p3 = midPoint;
	q2->p0 = midPoint;

	}
	*/
	return quads;
}


CVSystem::CD_Cdt CVSystem::Triangulator1::GetCDT(std::vector<CVSystem::MyPoint> pPoly, std::vector<std::vector<CVSystem::MyPoint>> cPoly, bool shouldRefine, bool isST)
{
	CD_Cdt cd_cdt;
	std::vector<std::vector<CD_VHandle>> handles;

	std::vector<CD_VHandle> hdls1;

	double lastX, lastY;
	lastX = lastY = -99;

	int pPolySize;
	if (isST)
	{
		pPolySize = pPoly.size();
	}
	else
	{
		pPolySize = pPoly.size() - 1;
	}
	for (size_t a = 0; a < pPolySize && pPolySize > 2; a++)
	{
		if (a == 0) {
			hdls1.push_back(cd_cdt.insert(CD_Point(pPoly[a].x, pPoly[a].y)));
		}
		else {
			if (!UtilityFunctions::pointEqual(lastX, lastY, pPoly[a].x, pPoly[a].y, 1e-5))
				hdls1.push_back(cd_cdt.insert(CD_Point(pPoly[a].x, pPoly[a].y)));
		}
		lastX = pPoly[a].x;
		lastY = pPoly[a].y;
	}
	handles.push_back(hdls1);

	if (cPoly.size() > 0)
	{
		for (size_t a = 0; a < cPoly.size(); a++)
		{
			std::vector<CVSystem::MyPoint> cPs = cPoly[a];
			std::vector<CD_VHandle> hdls2;

			lastX, lastY;
			lastX = lastY = -99;

			int cPsSize;
			if (isST)
			{
				cPsSize = cPs.size();
			}
			else
			{
				cPsSize = cPs.size() - 1;
			}
			for (size_t b = 0; b < cPsSize & cPsSize > 2; b++)
			{
				if (b == 0) {
					hdls2.push_back(cd_cdt.insert(CD_Point(cPs[b].x, cPs[b].y)));
				}
				else {
					if (!UtilityFunctions::pointEqual(lastX, lastY, cPs[b].x, cPs[b].y, 1e-5))
						hdls2.push_back(cd_cdt.insert(CD_Point(cPs[b].x, cPs[b].y)));
				}
				lastX = cPs[b].x;
				lastY = cPs[b].y;
			}
			handles.push_back(hdls2);
		}
	}
	// Get holes
	std::list<CD_Point> seeds = GetHoleSeeds(pPoly, cPoly, isST);
	// set constraints
	for (size_t a = 0; a < handles.size(); a++)
	{
		std::vector<CD_VHandle> hdls = handles[a];
		for (size_t b = 0; b < hdls.size(); b++)
		{
			if (b == 0) { cd_cdt.insert_constraint(hdls[hdls.size() - 1], hdls[b]); }
			else { cd_cdt.insert_constraint(hdls[b - 1], hdls[b]); }
		}
	}
	if (shouldRefine) { CGAL::refine_Delaunay_mesh_2(cd_cdt, seeds.begin(), seeds.end(), CD_Criteria(SystemParams::t_delaunay_aspect_bound, SystemParams::t_delaunay_max_length)); }
	else
	{	// use small factor to obtain less Delaunay
		CGAL::refine_Delaunay_mesh_2(cd_cdt, seeds.begin(), seeds.end(), CD_Criteria(0.01, SystemParams::t_delaunay_max_length));
	}
	return cd_cdt;
}

CVSystem::CD_Cdt CVSystem::Triangulator1::GetSTCDTWithRectangleBorder(std::vector<std::vector<MyPoint>> polygons, std::vector<std::vector<MyPoint>> STpolygons)
{
	CD_Cdt cd_cdt;
	std::vector<std::vector<CD_VHandle>> vHandles;

	std::vector<MyPoint> border;
	border.push_back(MyPoint(0, 0));
	border.push_back(MyPoint(this->_w_scaled, 0));
	border.push_back(MyPoint(this->_w_scaled, this->_h_scaled));
	border.push_back(MyPoint(0, this->_h_scaled));

	//polygons.push_back(border);
	for (size_t a = 0; a < polygons.size(); a++)
	{
		std::vector<CVSystem::MyPoint> cPs = polygons[a];
		std::vector<CD_VHandle> hdls;
		double lastX, lastY;
		lastX = lastY = -99;
		for (size_t b = 0; b < cPs.size() - 1 && cPs.size() - 1 > 2; b++)
		{
			if (b == 0) {
				hdls.push_back(cd_cdt.insert(CD_Point(cPs[b].x, cPs[b].y)));
			}
			else {
				if (!UtilityFunctions::pointEqual(lastX, lastY, cPs[b].x, cPs[b].y, 1e-5))
					hdls.push_back(cd_cdt.insert(CD_Point(cPs[b].x, cPs[b].y)));
			}
			lastX = cPs[b].x;
			lastY = cPs[b].y;
		}
		vHandles.push_back(hdls);
	}

	for (size_t a = 0; a < STpolygons.size(); a++)
	{
		std::vector<CVSystem::MyPoint> cPs = STpolygons[a];
		std::vector<CD_VHandle> hdls;
		double lastX, lastY;
		lastX = lastY = -99;
		for (size_t b = 0; b < cPs.size() && cPs.size() > 2; b++)
		{
			if (b == 0) {
				hdls.push_back(cd_cdt.insert(CD_Point(cPs[b].x, cPs[b].y)));
			}
			else {
				if (!UtilityFunctions::pointEqual(lastX, lastY, cPs[b].x, cPs[b].y, 1e-5))
					hdls.push_back(cd_cdt.insert(CD_Point(cPs[b].x, cPs[b].y)));
			}
			lastX = cPs[b].x;
			lastY = cPs[b].y;
		}
		//vHandles.push_back(hdls);
	}
	std::cout << "FK\n";
	// For image border only
	std::vector<CD_VHandle> borderVHandles;
	for (size_t a = 0; a < border.size(); a++) { borderVHandles.push_back(cd_cdt.insert(CD_Point(border[a].x, border[a].y))); }
	vHandles.push_back(borderVHandles);
	for (size_t a = 0; a < vHandles.size(); a++)
	{
		std::vector<CD_VHandle> hdls = vHandles[a];
		for (size_t b = 0; b < hdls.size(); b++)
		{
			if (b == 0) { cd_cdt.insert_constraint(hdls[hdls.size() - 1], hdls[b]); }
			else { cd_cdt.insert_constraint(hdls[b - 1], hdls[b]); }
		}
	}
	CGAL::refine_Delaunay_mesh_2(cd_cdt, CD_Criteria(SystemParams::t_delaunay_aspect_bound, SystemParams::t_delaunay_max_length));
	return cd_cdt;
}



// Add index on vertices similar to .OBJ file format
void CVSystem::Triangulator1::AssignVertexIndex(CVSystem::CD_Cdt& cdt)
{
	// Make all invalid
	for (CD_Cdt::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
	{
		vit->info() = -1;
	}//{ vit->info() = -2; }

	for (CD_Cdt::Finite_faces_iterator fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit)
	{
		if (fit->is_in_domain()) { fit->vertex(0)->info() = -1; fit->vertex(1)->info() = -1; fit->vertex(2)->info() = -1; }
	}

	// Start assigning index
	int idxCounter = 0;
	for (CD_Cdt::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
	{
		if (vit->info() == -1) { vit->info() = idxCounter++; }
	}

	// Add indexed vertices to list
	this->_vertexList.clear();
	for (CD_Cdt::Finite_vertices_iterator vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit)
		this->_vertexList.push_back(MyPoint(vit->point().x(), vit->point().y()));
	//{ if(vit->info() >= 0) { this->_vertexList.push_back(MyPoint(vit->point().x(), vit->point().y())); } }
}

std::vector<CVSystem::MyPoint> CVSystem::Triangulator1::GetPolygonLS(potrace_path_t* path, bool isDraw)
{
	potrace_curve_t curve = path->curve;
	std::vector<MyPoint> poly;

	for (int a = 0; a < curve.n; a++)
	{
		MyPoint p0; MyPoint p1; MyPoint p2; MyPoint p3;

		if (a == 0) { p0.x = curve.c[curve.n - 1][2].x; p0.y = this->_h_scaled - curve.c[curve.n - 1][2].y; }
		else { p0.x = curve.c[a - 1][2].x;       p0.y = this->_h_scaled - curve.c[a - 1][2].y; }

		p3 = MyPoint(curve.c[a][2].x, this->_h_scaled - curve.c[a][2].y);

		if (curve.tag[a] == POTRACE_CORNER)
		{
			MyPoint pMid(curve.c[a][1].x, this->_h_scaled - curve.c[a][1].y);

			CurveInterpolation::PointInterpolation(poly, p0, pMid, 0.5, SystemParams::t_subdivide_limit);
			poly.push_back(pMid);
			CurveInterpolation::PointInterpolation(poly, pMid, p3, 0.5, SystemParams::t_subdivide_limit);
		}
		else
		{
			p1 = MyPoint(curve.c[a][0].x, this->_h_scaled - curve.c[a][0].y);
			p2 = MyPoint(curve.c[a][1].x, this->_h_scaled - curve.c[a][1].y);

			CurveInterpolation::DeCasteljau(poly, p0, p1, p2, p3, SystemParams::t_subdivide_limit);
		}

		if (a == curve.n - 1) { poly.push_back(p3); }
	}
	if (isDraw)
	{
		for (int i = 0; i < poly.size(); i++)
		{
			if (path->sign == '+')
				cv::circle(traceImg, cv::Point2d(poly[i].x, poly[i].y), 1, cv::Scalar(255, 0, 255), -1);
			else
				cv::circle(traceImg, cv::Point2d(poly[i].x, poly[i].y), 1, cv::Scalar(0, 255, 0), -1);
		}
	}
	///***
	if (this->_has_tone) { _cssSmoothing->SmoothCurve2(poly, this->_distance_img_scaled); } // 
	if (isDraw)
	{
		for (int i = 0; i < poly.size(); i++)
		{
			if (path->sign == '+')
				cv::circle(traceImgS, cv::Point2d(poly[i].x, poly[i].y), 1, cv::Scalar(255, 0, 255), -1);
			else
				cv::circle(traceImgS, cv::Point2d(poly[i].x, poly[i].y), 1, cv::Scalar(0, 255, 0), -1);
		}
	}
	CleanUpPoly(poly);

	// !!! bug here...
	if (poly.size() > 2)
	{
		double dist = poly[0].Distance(poly[poly.size() - 1]);
		if (dist > M_EPS) { poly.push_back(poly[0]); }
	}

	return poly;
}

void CVSystem::Triangulator1::CurveToBezier(MyPoint p0, MyPoint p1, MyPoint p2, MyPoint p3, MyPoint& cp0, MyPoint& cp1)
{
	double xc1 = (p0.x + p1.x) / 2.0;		double yc1 = (p0.y + p1.y) / 2.0;
	double xc2 = (p1.x + p2.x) / 2.0;		double yc2 = (p1.y + p2.y) / 2.0;
	double xc3 = (p2.x + p3.x) / 2.0;		double yc3 = (p2.y + p3.y) / 2.0;

	double len1 = sqrt((p1.x - p0.x) * (p1.x - p0.x) + (p1.y - p0.y) * (p1.y - p0.y));
	double len2 = sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
	double len3 = sqrt((p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y));

	double k1 = len1 / (len1 + len2);		double k2 = len2 / (len2 + len3);

	double xm1 = xc1 + (xc2 - xc1) * k1;	double ym1 = yc1 + (yc2 - yc1) * k1;
	double xm2 = xc2 + (xc3 - xc2) * k2;	double ym2 = yc2 + (yc3 - yc2) * k2;

	// Resulting control points. Here smooth_value is mentioned
	// above coefficient K whose value should be in range [0...1].
	cp0.x = xm1 + (xc2 - xm1) * SystemParams::t_smooth_factor + p1.x - xm1;
	cp0.y = ym1 + (yc2 - ym1) * SystemParams::t_smooth_factor + p1.y - ym1;

	cp1.x = xm2 + (xc2 - xm2) * SystemParams::t_smooth_factor + p2.x - xm2;
	cp1.y = ym2 + (yc2 - ym2) * SystemParams::t_smooth_factor + p2.y - ym2;
}

potrace_bitmap_t* CVSystem::Triangulator1::GetBM(cv::Mat img)
{
	cv::Mat flipImg;
	cv::flip(img, flipImg, 0);
	cv::Mat inpImg;
	flipImg.convertTo(inpImg, CV_8UC1);

	potrace_bitmap_t* bm = bm_new(img.cols, img.rows);
	for (int a = 0; a < bm->w; a++)
	{
		for (int b = 0; b < bm->h; b++)
		{
			// Todo: users should be able to set the threshold
			if (inpImg.ptr<uchar>(b, a)[0] < 127) { BM_PUT(bm, a, b, 1); }
			else { BM_PUT(bm, a, b, 0); }
		}
	}
	return bm;
}

// LS
std::vector<CVSystem::MyQuad> CVSystem::Triangulator1::GetLSQuads()
{
	double val = 1.0 / SystemParams::t_scale_factor;
	std::vector<CVSystem::MyQuad> resizedLSQuads;
	for (size_t a = 0; a < _lsQuads.size(); a++)
	{
		resizedLSQuads.push_back(_lsQuads[a].Resize(val));
	}
	return resizedLSQuads;
}

std::vector<CVSystem::MyPoint> CVSystem::Triangulator1::GetVertexList()
{
	double val = 1.0 / SystemParams::t_scale_factor;
	std::vector<CVSystem::MyPoint>	resizedVertexList;
	for (size_t a = 0; a < this->_vertexList.size(); a++)
	{
		resizedVertexList.push_back(this->_vertexList[a].Resize(val));
	}
	return resizedVertexList;
}

std::vector<CVSystem::MyIndexedTriangle> CVSystem::Triangulator1::GetIndexedTriangles()
{
	return this->_indexedTriangles;
}

std::vector<CVSystem::MyIndexedTriangle> CVSystem::Triangulator1::GetBorderWTriangles()
{
	return this->_borderWTriangles;
}

std::vector<CVSystem::MyIndexedTriangle> CVSystem::Triangulator1::GetBorderSCTriangles()
{
	return this->_borderSCTriangles;
}

std::vector<CVSystem::MyIndexedTriangle> CVSystem::Triangulator1::GetWTriangles()
{
	return this->_wTriangles;
}

std::vector<CVSystem::MyIndexedTriangle> CVSystem::Triangulator1::GetSCTriangles()
{
	return this->_scTriangles;
}

std::vector<CVSystem::MyIndexedBezierCurves> CVSystem::Triangulator1::GetIndexedBezierCurves()
{
	for (size_t a = 0; a < _indexedBezierCurves.size(); a++)
	{
		_indexedBezierCurves[a].CalculateOrientation(_vertexList);
	}

	return this->_indexedBezierCurves;
}

std::vector<CVSystem::MyIndexedTriangle> CVSystem::Triangulator1::GetBTriangles()
{
	return this->_bTriangles;
}

// Deprecated
// Get corner points
std::vector<CVSystem::MyPoint> CVSystem::Triangulator1::GetCornerList()
{
	double val = 1.0 / SystemParams::t_scale_factor;
	std::vector<CVSystem::MyPoint>	resizedCornerList;
	for (size_t a = 0; a < this->_cornerList.size(); a++)
	{
		resizedCornerList.push_back(this->_cornerList[a].Resize(val));
	}
	return resizedCornerList;
}

std::vector<int> CVSystem::Triangulator1::GetCornerIndices()
{
	return this->_cornerIndices;
}

std::vector<CVSystem::MyIndexedLine> CVSystem::Triangulator1::GetCornerEdges()
{
	return this->_cornerEdges;
}
// Deprecated
std::vector<CVSystem::MyPoint> CVSystem::Triangulator1::GetPolygonO(potrace_curve_t curve)
{
	std::vector<CVSystem::MyPoint> pts;
	for (int a = 0; a < curve.n; a++)
	{
		if (a == 0) { pts.push_back(CVSystem::MyPoint(curve.c[curve.n - 1][2].x, this->_h_scaled - curve.c[curve.n - 1][2].y)); }
		else { pts.push_back(CVSystem::MyPoint(curve.c[a - 1][2].x, this->_h_scaled - curve.c[a - 1][2].y)); }

		if (curve.tag[a] == POTRACE_CORNER) { pts.push_back(CVSystem::MyPoint(curve.c[a][1].x, this->_h_scaled - curve.c[a][1].y)); }
	}
	return pts;
}

// Deprecated
std::vector<CVSystem::MyQuad> CVSystem::Triangulator1::GetBezierQuad0(potrace_path_t* path)
{
	std::vector<CVSystem::MyQuad> quads;
	potrace_curve_t curve = path->curve;
	bool isPlus = (path->sign == '-') ? false : true;

	// build polygon
	std::vector<CD_Point> points;
	for (int a = 0; a < curve.n; a++)
	{
		CD_Point p(curve.c[a][2].x, this->_h_scaled - curve.c[a][2].y);
		points.push_back(p);
	}

	for (int a = 0; a < curve.n; a++)
	{
		if (curve.tag[a] == POTRACE_CORNER) {}
		else // curve
		{
			MyPoint p0; MyPoint p1; MyPoint p2; MyPoint p3;

			MyPoint pPrev; MyPoint pNext;

			if (a == 0)
			{
				p0 = MyPoint(curve.c[curve.n - 1][2].x, this->_h_scaled - curve.c[curve.n - 1][2].y);
				pPrev = MyPoint(curve.c[curve.n - 2][2].x, this->_h_scaled - curve.c[curve.n - 2][2].y);
			}
			else
			{
				p0 = MyPoint(curve.c[a - 1][2].x, this->_h_scaled - curve.c[a - 1][2].y);
				if (a == 1)
				{
					pPrev = MyPoint(curve.c[curve.n - 1][2].x, this->_h_scaled - curve.c[curve.n - 1][2].y);
				}
				else
				{
					pPrev = MyPoint(curve.c[a - 2][2].x, this->_h_scaled - curve.c[a - 2][2].y);
				}
			}

			p3 = MyPoint(curve.c[a][2].x, this->_h_scaled - curve.c[a][2].y);
			if (a == curve.n - 1)
			{
				pNext = MyPoint(curve.c[0][2].x, this->_h_scaled - curve.c[0][2].y);
			}
			else
			{
				pNext = MyPoint(curve.c[a + 1][2].x, this->_h_scaled - curve.c[a + 1][2].y);
			}

			CurveToBezier(pPrev, p0, p3, pNext, p1, p2);

			// Using data from Potrace
			/*p1.x = curve.c[a][0].x;
			p1.y = h - curve.c[a][0].y;

			p2.x = curve.c[a][1].x;
			p2.y = h - curve.c[a][1].y;*/

			quads.push_back(CVSystem::MyQuad(p0, p1, p2, p3));
		}
	}
	return quads;
}
std::list<CVSystem::CD_Point> CVSystem::Triangulator1::GetHoleSeeds(std::vector<CVSystem::MyPoint> pPoly, std::vector<std::vector<CVSystem::MyPoint>> cPoly, bool isST)
{
	std::list<CD_Point> seeds;

	// (PD) for outer polygon
	PD_Cdt pd_cdt;
	Polygon_2 polygon1;

	int pPolySize;
	if (isST)
	{
		pPolySize = pPoly.size();
	}
	else
	{
		pPolySize = pPoly.size() - 1;
	}
	for (size_t a = 0; a < pPolySize; a++)
	{
		polygon1.push_back(PD_Point(pPoly[a].x, pPoly[a].y));
	}
	insert_polygon(pd_cdt, polygon1);

	// (PD )for inner polygon
	for (size_t a = 0; a < cPoly.size(); a++)
	{
		Polygon_2 polygon2;
		std::vector<CVSystem::MyPoint> cPs = cPoly[a];

		int cPsSize;
		if (isST)
		{
			cPsSize = cPs.size();
		}
		else
		{
			cPsSize = cPs.size() - 1;
		}
		for (size_t b = 0; b < cPsSize; b++)
		{
			polygon2.push_back(PD_Point(cPs[b].x, cPs[b].y));
		}
		insert_polygon(pd_cdt, polygon2);
	}

	// (PD) Mark domains
	mark_domains(pd_cdt);

	// (PD) Get holes
	for (PD_Cdt::Finite_faces_iterator fit = pd_cdt.finite_faces_begin(); fit != pd_cdt.finite_faces_end(); ++fit)
	{
		if (!fit->info().in_domain())
		{
			float x_center = ((*(fit->vertex(0))).point().x() + (*(fit->vertex(1))).point().x() + (*(fit->vertex(2))).point().x()) / 3.0;
			float y_center = ((*(fit->vertex(0))).point().y() + (*(fit->vertex(1))).point().y() + (*(fit->vertex(2))).point().y()) / 3.0;
			seeds.push_back(CD_Point(x_center, y_center));
		}
	}
	return seeds;
}

void CVSystem::Triangulator1::Resubdivide(std::vector<CVSystem::MyPoint>& poly)
{
	using namespace CVSystem;
	std::vector<CVSystem::MyPoint> tempPoly;
	for (size_t a = 0; a < poly.size(); a++)
	{
		tempPoly.push_back(MyPoint(poly[a].x, poly[a].y));
	}

	poly.clear();

	// Warning, hard parameters, but work just fine right now...
	double factr = 0.25;
	if (SystemParams::t_subdivide_limit > 5.0) { factr = 0.01; }

	for (size_t a = 0; a < tempPoly.size() - 1; a++)
	{
		CurveInterpolation::PointInterpolation(poly, tempPoly[a], tempPoly[a + 1], 0.5, SystemParams::t_subdivide_limit * factr);
		poly.push_back(tempPoly[a + 1]);
	}
}

void CVSystem::Triangulator1::RDPSimplification(std::vector<CVSystem::MyPoint> poly, std::vector<CVSystem::MyPoint>& rdpPoly, std::vector<bool>& bFlag)
{
	bFlag = std::vector<bool>(poly.size());
	for (size_t a = 0; a < poly.size(); a++) { bFlag[a] = false; }
	bFlag[0] = true;
	bFlag[poly.size() - 1] = true;

	CurveRDP::RDP(bFlag, poly, SystemParams::t_rdp_epsilon, 0, poly.size() - 1, SystemParams::t_rdp_point_min);
	for (size_t a = 0; a < poly.size(); a++)
	{
		if (bFlag[a]) rdpPoly.push_back(poly[a]);
	}

	//CleanUpPoly(rdpPoly);
}


int CVSystem::Triangulator1::ReduceSTpoint(std::vector<std::vector<MyPoint>> polygons, std::vector<CVSystem::MyPoint> souPoly, std::vector<CVSystem::MyPoint>& desPoly)
{
	int reduceCount = 0;
	std::vector<MyPoint> tNewPoly;
	std::vector<MyPoint> newPoly;

	std::vector<int> tNewPolySource; //-1 st 0~n ori
	std::vector<int> tNewPolyIdx;

	std::vector<PairData> PolyIdxCount;

	double lastX, lastY;
	lastX = lastY = -99;


	int oriSource = -2;
	for (int p = 0; p < souPoly.size(); p++)
	{
		MyPoint tp = souPoly[p];
		MyPoint mindisP;
		int PolyIdx;
		int PolySimpleIdx;
		double dis = 99;
		for (int i = 0; i < polygons.size(); i++)
		{
			for (int j = 0; j < polygons[i].size() - 1; j++)
			{
				double tempDis = tp.Distance(polygons[i][j]);
				if (tempDis < dis)
				{
					mindisP = polygons[i][j];
					PolyIdx = i;
					PolySimpleIdx = j;
					dis = tempDis;
				}
			}
		}

		if (dis > SystemParams::t_STDis)
		{
			if (!UtilityFunctions::pointEqual(lastX, lastY, tp.x, tp.y, 1e-3))
			{
				bool findCC = false;
				for (int CC = 0; CC < PolyIdxCount.size(); CC++)
				{
					if (PolyIdxCount[CC].first == -1)
					{
						PolyIdxCount[CC].second++;
						findCC = true;
					}
				}
				if (!findCC)
					PolyIdxCount.push_back(PairData(-1, 1));
			}
		}
		else
		{
			if (!UtilityFunctions::pointEqual(lastX, lastY, mindisP.x, mindisP.y, 1e-3))
			{
				bool findCC = false;
				for (int CC = 0; CC < PolyIdxCount.size(); CC++)
				{
					if (PolyIdxCount[CC].first == PolyIdx)
					{
						PolyIdxCount[CC].second++;
						findCC = true;
					}
				}
				if (!findCC)
					PolyIdxCount.push_back(PairData(PolyIdx, 1));
			}

		}
	}

	UtilityFunctions::QuickSortMain(PolyIdxCount);
	if (PolyIdxCount[PolyIdxCount.size() - 1].first != -1)
		oriSource = PolyIdxCount[PolyIdxCount.size() - 1].first;
	else
		oriSource = PolyIdxCount[PolyIdxCount.size() - 2].first;

	for (int p = 0; p < souPoly.size(); p++)
	{
		MyPoint tp = souPoly[p];
		MyPoint mindisP;
		int PolyIdx;
		int PolySimpleIdx;
		double dis = 99;
		for (int j = 0; j < polygons[oriSource].size() - 1; j++)
		{
			double tempDis = tp.Distance(polygons[oriSource][j]);
			if (tempDis < dis)
			{
				mindisP = polygons[oriSource][j];
				PolyIdx = oriSource;
				PolySimpleIdx = j;
				dis = tempDis;
			}
		}

		if (dis > SystemParams::t_STDis)
		{
			if (!UtilityFunctions::pointEqual(lastX, lastY, tp.x, tp.y, 1e-3))
			{
				if ((tNewPoly.size() > 0 && !UtilityFunctions::pointEqual(tNewPoly[0].x, tNewPoly[0].y, tp.x, tp.y, 1e-3))
					|| tNewPoly.size() == 0)
				{
					tNewPoly.push_back(tp);
					tNewPolySource.push_back(-1);
					tNewPolyIdx.push_back(p);
					lastX = tp.x;
					lastY = tp.y;
				}
			}
		}
		else
		{
			if (!UtilityFunctions::pointEqual(lastX, lastY, mindisP.x, mindisP.y, 1e-3))
			{
				if ((tNewPoly.size() > 0 && !UtilityFunctions::pointEqual(tNewPoly[0].x, tNewPoly[0].y, mindisP.x, mindisP.y, 1e-3))
					|| tNewPoly.size() == 0)
				{
					tNewPoly.push_back(mindisP);
					tNewPolySource.push_back(PolyIdx);
					tNewPolyIdx.push_back(PolySimpleIdx);
					lastX = mindisP.x;
					lastY = mindisP.y;

					reduceCount++;
				}
			}

		}
	}


	std::cout << SystemParams::t_STDis << "  " << souPoly.size() << "  " << reduceCount << std::endl;

	cv::Mat cc = testImg.clone();
	std::stringstream ss;

	for (int Idx = 0; Idx < souPoly.size(); Idx++)
	{
		cv::circle(cc, cv::Point2d(souPoly[Idx].x, souPoly[Idx].y), 2, cv::Scalar(255, 0, 0), -1);
	}
	for (int Idx = 0; Idx < tNewPoly.size(); Idx++)
	{
		cv::circle(cc, cv::Point2d(tNewPoly[Idx].x, tNewPoly[Idx].y), 1, cv::Scalar(0, 255, 0), -1);
	}
	ss << "rdPoly_" << testcount << "_t" << 1 << "_" << SystemParams::t_STDis << ".png";
	cv::imwrite(ss.str().c_str(), cc);

	int ooriS = oriSource;
	testcount++;
	if (false && testcount != 9 && testcount != 11 && testcount != 48 && testcount != 50)
	{
		std::cout << testcount << std::endl;
		int ooriS = oriSource;
		if (souPoly.size() == reduceCount)
		{
			desPoly.clear();
			for (int j = 0; j < polygons[oriSource].size() - 1; j++)
			{
				desPoly.push_back(polygons[oriSource][j]);
			}
			return 0;
		}
		else
		{
			bool isW = false;//debug
			int step = 0;
			for (int p = 0; p < tNewPoly.size(); p++)
			{
				//分析曲線
				int oriPathCount = 0;
				for (int i = -2; i < 3; i++)
				{
					int pIdx = p + i;
					if (pIdx < 0)
					{
						pIdx += tNewPoly.size();
					}
					pIdx = (pIdx >= tNewPoly.size()) ? pIdx - tNewPoly.size() : pIdx;
					if (tNewPolySource[pIdx] != -1)
					{
						oriPathCount++;
						oriSource = tNewPolySource[pIdx];
					}
				}
				if (oriPathCount >= 3)//貼合原始曲線
				{
					int last = -1;
					int next = -1;
					int lastPP = -1;
					int nextPP = -1;
					for (int i = -2; i < 3; i++)
					{
						int pIdx = p + i;
						if (pIdx < 0)
						{
							pIdx += tNewPoly.size();
						}
						if (pIdx >= tNewPoly.size())
						{
							pIdx -= tNewPoly.size();
						}
						if (tNewPolySource[pIdx] != -1)
						{
							if (i < 0)
							{
								last = tNewPolyIdx[pIdx];
								lastPP = pIdx;
							}
							else if (i >= 0 && next == -1)
							{
								next = tNewPolyIdx[pIdx];
								nextPP = pIdx;
							}
						}
					}

					for (int eraseIdx = lastPP + 1; eraseIdx < nextPP; eraseIdx++)
					{
						if (tNewPolySource[eraseIdx] == -1)
						{
							tNewPoly.erase(tNewPoly.begin() + eraseIdx);
							tNewPolySource.erase(tNewPolySource.begin() + eraseIdx);
							tNewPolyIdx.erase(tNewPolyIdx.begin() + eraseIdx);
						}
					}

					if (last != -1 && abs(last - next) != 1)
					{
						if (last > next && (last - next) < polygons[oriSource].size() / 2)//new
						{
							tNewPolyIdx[lastPP] = next;
							tNewPolyIdx[nextPP] = last;

							MyPoint tp = tNewPoly[lastPP];
							tNewPoly[lastPP] = tNewPoly[nextPP];
							tNewPoly[nextPP] = tp;

							int t = last;
							last = next;
							next = t;
						}
						int insertIdx = last;
						insertIdx++;
						insertIdx = (insertIdx >= (polygons[oriSource].size() - 1)) ? insertIdx - (polygons[oriSource].size() - 1) : insertIdx;
						for (; insertIdx != last && insertIdx != next;)
						{
							MyPoint tp = polygons[oriSource][insertIdx];
							bool isExist = false;
							for (int findIdx = 0; findIdx < tNewPoly.size(); findIdx++)
							{
								if (UtilityFunctions::pointEqual(tNewPoly[findIdx].x, tNewPoly[findIdx].y, tp.x, tp.y, 1e-3))
								{
									isExist = true;
									break;
								}
							}
							if (!isExist)
							{
								tNewPoly.insert(tNewPoly.begin() + p, tp);
								tNewPolySource.insert(tNewPolySource.begin() + p, oriSource);
								tNewPolyIdx.insert(tNewPolyIdx.begin() + p, insertIdx);
								p++;

								/*
								cv::Mat ccc = testImg.clone();
								cv::circle(ccc, cv::Point2d(polygons[oriSource][last].x, polygons[oriSource][last].y), 2, cv::Scalar(0, 0, 100), -1);
								cv::circle(ccc, cv::Point2d(polygons[oriSource][next].x, polygons[oriSource][next].y), 2, cv::Scalar(0, 0, 200), -1);
								cv::circle(ccc, cv::Point2d(tNewPoly[p-1].x, tNewPoly[p-1].y), 2, cv::Scalar(0, 255, 255), -1);
								for(int i = -2; i < 3; i++)
								{
								int R = (p + i < 0) ? tNewPoly.size() + p + i : p + i;
								R = (R >= tNewPoly.size()) ? R - tNewPoly.size() : R;

								cv::circle(ccc, cv::Point2d(tNewPoly[R].x, tNewPoly[R].y), 1, cv::Scalar(0, 255, 0), -1);
								ss.str("");
								ss<<"rdPoly_"<<testcount<<"_t"<<2<<"_"<<SystemParams::t_STDis<<"_"<<p-1<<"_"<<R<<"_"<<lastPP<<"_"<<nextPP<<"_"<<last<<"_"<<next<<".png";

								cv::imwrite(ss.str().c_str(), ccc);
								}
								*/
							}
							insertIdx++;
							insertIdx = (insertIdx >= (polygons[oriSource].size() - 1)) ? insertIdx - (polygons[oriSource].size() - 1) : insertIdx;
						}
					}

				}
				else {
					if (tNewPolySource[p] != -1)
					{
						int last = -1;
						int next = -1;
						int lastPP = -1;
						int nextPP = -1;
						for (int i = -2; i < 3; i++)
						{
							int pIdx = p + i;
							if (pIdx < 0)
							{
								pIdx = tNewPoly.size() + pIdx;
							}
							pIdx = (pIdx >= tNewPoly.size()) ? pIdx - tNewPoly.size() : pIdx;
							if (tNewPolySource[pIdx] != -1)
							{
								if (i < 0)
								{
									last = tNewPolyIdx[pIdx];
									lastPP = pIdx;
								}
								else if (i >= 0 && next == -1)
								{
									next = tNewPolyIdx[pIdx];
									nextPP = pIdx;
								}
							}
						}

						for (int eraseIdx = lastPP + 1; eraseIdx < nextPP; eraseIdx++)
						{
							if (tNewPolySource[eraseIdx] != -1)
							{
								tNewPoly.erase(tNewPoly.begin() + eraseIdx);
								tNewPolySource.erase(tNewPolySource.begin() + eraseIdx);
								tNewPolyIdx.erase(tNewPolyIdx.begin() + eraseIdx);
							}
						}

						if (last != -1)
						{
							int insertIdx = last;
							insertIdx++;
							insertIdx = (insertIdx >= souPoly.size()) ? insertIdx - souPoly.size() : insertIdx;
							for (; insertIdx != last && insertIdx != next;)
							{
								MyPoint tp = souPoly[insertIdx];
								bool isExist = false;
								for (int findIdx = 0; findIdx < tNewPoly.size(); findIdx++)
								{
									if (UtilityFunctions::pointEqual(tNewPoly[findIdx].x, tNewPoly[findIdx].y, tp.x, tp.y, 1e-3))
									{
										isExist = true;
										break;
									}
								}
								if (!isExist)
								{
									isW = true;
									cv::circle(cc, cv::Point2d(tp.x, tp.y), 1, cv::Scalar(255, 255, 0), -1);

									tNewPoly.insert(tNewPoly.begin() + p, tp);
									tNewPolySource.insert(tNewPolySource.begin() + p, -1);
									tNewPolyIdx.insert(tNewPolyIdx.begin() + p, insertIdx);
									p++;
								}

								insertIdx++;
								insertIdx = (insertIdx >= souPoly.size()) ? insertIdx - souPoly.size() : insertIdx;
							}
						}
					}
				}
			}
		}
	}
	/*
	cc = testImg.clone();
	for(int Idx = 0; Idx < souPoly.size(); Idx++)
	{
	cv::circle(cc, cv::Point2d(souPoly[Idx].x, souPoly[Idx].y), 2, cv::Scalar(255, 0, 0), -1);
	}
	bool hasST = false;
	for(int Idx = 0; Idx < tNewPoly.size(); Idx++)
	{
	cv::circle(cc, cv::Point2d(tNewPoly[Idx].x, tNewPoly[Idx].y), 1, cv::Scalar(0, 255, 0), -1);
	cv::Mat ccc = testImg.clone();
	cv::circle(ccc, cv::Point2d(tNewPoly[Idx].x, tNewPoly[Idx].y), 2, cv::Scalar(0, 255, 0), -1);
	ss.str("");
	ss<<"rdPoly_"<<testcount<<"_t"<<2<<"_"<<SystemParams::t_STDis<<"_"<<Idx<<"_"<<tNewPolySource[Idx]<<"_"<<tNewPolyIdx[Idx]<<".png";
	//cv::imwrite(ss.str().c_str(), ccc);

	if(tNewPolySource[Idx] == -1)
	hasST = true;
	}
	*/
	std::cout << "FK" << std::endl;
	bool hasST = false;
	for (int Idx = 0; Idx < tNewPoly.size(); Idx++)
	{
		if (tNewPolySource[Idx] == -1)
			hasST = true;
	}

	desPoly.clear();
	if (hasST)
	{



		cc = testImg.clone();
		for (int Idx = 0; Idx < souPoly.size(); Idx++)
		{
			cv::circle(cc, cv::Point2d(souPoly[Idx].x, souPoly[Idx].y), 2, cv::Scalar(255, 0, 0), -1);
		}
		for (int Idx = 0; Idx < tNewPoly.size(); Idx++)
		{
			cv::circle(cc, cv::Point2d(tNewPoly[Idx].x, tNewPoly[Idx].y), 1, cv::Scalar(0, 255, 0), -1);
			cv::Mat ccc = testImg.clone();
			cv::circle(ccc, cv::Point2d(tNewPoly[Idx].x, tNewPoly[Idx].y), 2, cv::Scalar(0, 255, 0), -1);
			ss.str("");
			ss << "rdPoly_" << testcount << "_t" << 2 << "_" << SystemParams::t_STDis << "_" << Idx << "_" << tNewPolySource[Idx] << "_" << tNewPolyIdx[Idx] << ".png";

			//cv::imwrite(ss.str().c_str(), ccc);

		}


		ss.str("");
		ss << "rdPoly_" << testcount++ << "_t" << 2 << "_" << SystemParams::t_STDis << ".png";
		cv::imwrite(ss.str().c_str(), cc);


		desPoly = tNewPoly;
		/*
		for(int j = 0; j < polygons[ooriS].size()-1; j++)
		{
		desPoly.push_back(polygons[ooriS][j]);
		}
		*/
		return 1;
	}
	else
	{
		for (int j = 0; j < polygons[ooriS].size() - 1; j++)
		{
			desPoly.push_back(polygons[ooriS][j]);
		}
		return 0;
	}
}

void  CVSystem::Triangulator1::ProcessEdgesFromCDT(CVSystem::CD_Cdt cd_cdt)
{
	this->_borderSCTriangles.clear();
	this->_borderWTriangles.clear();

	std::vector<MyIndexedLine> unsortedLines;

	for (CD_Cdt::Finite_edges_iterator eit = cd_cdt.finite_edges_begin(); eit != cd_cdt.finite_edges_end(); ++eit)
	{
		const CD_Cdt::Face_handle& fh1 = eit->first;					// left?
		const CD_Cdt::Face_handle& fh2 = fh1->neighbor(eit->second);	// right?
		if (fh1->triangleType == TRIANGLE_BLACK && fh2->triangleType != TRIANGLE_BLACK)
		{
			int i = eit->second;
			int index0 = fh1->vertex(fh1->cw(i))->info();
			int index1 = fh1->vertex(fh1->ccw(i))->info();

			unsortedLines.push_back(MyIndexedLine(index0, index1));
		}
		// should be flipped so black triangles are in left side
		if (fh1->triangleType != TRIANGLE_BLACK && fh2->triangleType == TRIANGLE_BLACK)
		{
			int i = eit->second;
			int index0 = fh1->vertex(fh1->cw(i))->info();
			int index1 = fh1->vertex(fh1->ccw(i))->info();

			unsortedLines.push_back(MyIndexedLine(index1, index0));
			//unsortedLines.push_back(MyLine(x2, y2, x1, y1, index1, index0));
		}
		if (fh1->triangleType == TRIANGLE_BLACK && fh2->triangleType == TRIANGLE_WHITE)
		{
			this->_borderWTriangles.push_back(MyIndexedTriangle(fh1->vertex(0)->info(), fh1->vertex(1)->info(), fh1->vertex(2)->info(), TRIANGLE_WHITE));
		}
		else if (fh1->triangleType == TRIANGLE_WHITE && fh2->triangleType == TRIANGLE_BLACK)
		{
			this->_borderWTriangles.push_back(MyIndexedTriangle(fh2->vertex(0)->info(), fh2->vertex(1)->info(), fh2->vertex(2)->info(), TRIANGLE_WHITE));
		}
		else if (fh1->triangleType == TRIANGLE_BLACK && fh2->triangleType == TRIANGLE_SCREENTONE)
		{
			this->_borderSCTriangles.push_back(MyIndexedTriangle(fh1->vertex(0)->info(), fh1->vertex(1)->info(), fh1->vertex(2)->info(), TRIANGLE_SCREENTONE));
		}
		else if (fh1->triangleType == TRIANGLE_SCREENTONE && fh2->triangleType == TRIANGLE_BLACK)
		{
			this->_borderSCTriangles.push_back(MyIndexedTriangle(fh2->vertex(0)->info(), fh2->vertex(1)->info(), fh2->vertex(2)->info(), TRIANGLE_SCREENTONE));
		}
	}
	this->_indexedBezierCurves = LinesSorter::BuildBezierCurves2(unsortedLines, this->_vertexList);
}

void CVSystem::Triangulator1::GetCornersFromPolygon()
{
	double cornerAngle = 2.792526803; // 160
	for (size_t a = 0; a < this->_indexedBezierCurves.size(); a++)
	{
		MyIndexedBezierCurves bCurve = this->_indexedBezierCurves[a];
		//this->_cornerList.push_back(poly.points[0]);		
		for (size_t b = 0; b < bCurve.indices.size() - 1; b++)
		{
			MyPoint p0;
			if (b == 0) { p0 = _vertexList[bCurve.indices[bCurve.indices.size() - 2]]; }
			else { p0 = _vertexList[bCurve.indices[b - 1]]; }

			MyPoint p1 = _vertexList[bCurve.indices[b]];
			MyPoint p2 = _vertexList[bCurve.indices[b + 1]];
			MyPoint dir0 = p0 - p1;
			MyPoint dir1 = p2 - p1;

			//cos(angle) = dot_product / (a.len * b.len) 
			double angle = acos(dir0.Dot(dir1) / (dir0.Length() * dir1.Length()));

			if (angle < cornerAngle)
			{
				this->_cornerList.push_back(p1);	// Fucking delete this
				this->_cornerIndices.push_back(bCurve.indices[b]);
			}
		}
	}
}

void CVSystem::Triangulator1::ProcessPolygonCorners()
{
	// all corner flags are false
	for (size_t a = 0; a < this->_indexedBezierCurves.size(); a++)
	{
		MyIndexedBezierCurves* bCurve = &this->_indexedBezierCurves[a];
		bCurve->cornerFlags.resize(bCurve->indices.size());
		for (size_t b = 0; b < bCurve->cornerFlags.size(); b++)
		{
			bCurve->cornerFlags[b] = false;
		}
	}

	// Todo: put this in SystemParams
	double cornerAngle = 2.792526803; // 160
	for (size_t a = 0; a < this->_indexedBezierCurves.size(); a++)
	{
		MyIndexedBezierCurves* bCurve = &this->_indexedBezierCurves[a];
		for (size_t b = 0; b < bCurve->indices.size() - 1; b++)
		{
			MyPoint p0;
			if (b == 0) { p0 = _vertexList[bCurve->indices[bCurve->indices.size() - 2]]; }
			else { p0 = _vertexList[bCurve->indices[b - 1]]; }

			MyPoint p1 = _vertexList[bCurve->indices[b]];
			MyPoint p2 = _vertexList[bCurve->indices[b + 1]];
			MyPoint dir0 = p0 - p1;
			MyPoint dir1 = p2 - p1;

			double angle = acos(dir0.Dot(dir1) / (dir0.Length() * dir1.Length()));

			if (angle < cornerAngle) { bCurve->cornerFlags[b] = true; }
		}
	}

	for (size_t a = 0; a < this->_indexedBezierCurves.size(); a++)
	{
		MyIndexedBezierCurves bCurve = this->_indexedBezierCurves[a];
		for (size_t b = 0; b < bCurve.indices.size() - 1; b++)
		{
			if (bCurve.cornerFlags[b])
			{
				int prevIndex = -1;
				int nextIndex = -1;
				if (b == 0) { prevIndex = bCurve.indices.size() - 2; }
				else { prevIndex = b - 1; }

				if (b == bCurve.indices.size() - 2) { nextIndex = 0; }
				else { nextIndex = b + 1; }

				// always take previous
				_cornerEdges.push_back(MyIndexedLine(bCurve.indices[prevIndex], bCurve.indices[b]));

				// if next is also corner, don't take that
				if (!bCurve.cornerFlags[nextIndex])
				{
					_cornerEdges.push_back(MyIndexedLine(bCurve.indices[b], bCurve.indices[nextIndex]));
				}

			}
		}
	}
}

void CVSystem::Triangulator1::DiscardVertexOffset(float offsetX, float offsetY)
{
	for (size_t i = 0; i < _vertexList.size(); i++)
	{
		_vertexList[i].x += offsetX;
		_vertexList[i].y += offsetY;
	}
}

void CVSystem::Triangulator1::CleanUpPoly(std::vector<CVSystem::MyPoint>& poly)
{
	std::vector<CVSystem::MyPoint> tempPoly;
	for (size_t a = 0; a < poly.size(); a++)
	{
		tempPoly.push_back(CVSystem::MyPoint(poly[a].x, poly[a].y));
	}

	poly.clear();
	poly.push_back(tempPoly[0]);
	for (size_t a = 1; a < tempPoly.size(); a++)
	{
		if (tempPoly[a] != tempPoly[a - 1])
		{
			poly.push_back(tempPoly[a]);
		}
	}

	// get size
	std::vector<cv::Point> contour;
	for (size_t b = 0; b < poly.size(); b++)
	{
		contour.push_back(cv::Point2d(poly[b].x, poly[b].y));
	}
	double polyArea = cv::contourArea(contour);

	if (polyArea < SystemParams::t_min_poly_area * SystemParams::t_scale_factor)
	{ /*std::cout << "poly_cleared ";*/ poly.clear();
	}
}
