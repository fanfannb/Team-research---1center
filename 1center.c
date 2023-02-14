#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define GNUPLOT_PATH "C:\\PROGRA~1\\gnuplot\\bin\\gnuplot.exe"// gnuplot.exe的位置
#define RANGE "set xrange[-150:150]\nset yrange[-150:150]\n"// 建立坐标系大小
#define NMAX 100


struct triangle { int i, j, k; };
double x[NMAX], y[NMAX], z[NMAX];
struct triangle out[NMAX];
struct point { double x, y; };
struct point p[1000];
struct point arr[1000];
struct point vp[1000]; //voronoi 点数
struct point inst[100];//关于1中心问题的交叉点数

struct point intersection(struct point p1, struct point p2, struct point p3, struct point p4) {
	double m1, c1, m2, c2;
	double x1, y1, x2, y2;
	double dx, dy;
	struct point ins;


	dx = p2.x - p1.x;

	dy = p2.y - p1.y;

	m1 = dy / dx;
	// y = mx + c

	// 截距 c = y - mx

	c1 = p1.y - m1 * p1.x; // 如同   y2 - 斜率 * x2  一样



	dx = p4.x - p3.x;

	dy = p4.y - p3.y;

	m2 = dy / dx;

	// y = mx + c

	// 截距 c = y - mx

	c2 = p3.y - m2 * p3.x; // 如同  y2 - 斜率 * x2  一样
	
	if ((m1 - m2) == 0) {
		printf("这些线段之间没有交叉\n");
		ins.x = -999;
		ins.y = -999;
		return ins;
	}
	else if (p1.x - p2.x == 0) {

		ins.x = p1.x;
		ins.y = m2 * p1.x + c2;

		return ins;


	}
	else if (p3.x - p4.x == 0) {

		ins.x = p3.x;
		ins.y = m1 * p3.x + c1;

		return ins;


	}


	else {

		ins.x = (c2 - c1) / (m1 - m2);

		ins.y = m1 * ins.x + c1;

		//printf("交叉点: = %.2f, %.2f\n", intersection_X, intersection_Y);
		return ins;

	}


}


double area(struct point q1, struct point q2, struct point q3)
{
	return((q1.x - q3.x)*(q2.y - q3.y) + (q2.x - q3.x)*(q3.y - q1.y));
}


void quicksort(int first, int last)
{
	int i, j, somewhere;
	struct point pivot, temp;

	if (first < last)
	{
		somewhere = (first + last) / 2;
		pivot = p[somewhere];
		i = first; j = last;
		do {
			while (area(p[0], p[i], pivot) > 0)i++;
			while (area(p[0], p[j], pivot) < 0)j--;
			if (i <= j)
			{
				temp = p[i]; p[i] = p[j]; p[j] = temp;
				i++; j--;
			}
		} while (i <= j);

		quicksort(first, j);
		quicksort(i, last);
	}
}
void quicksort2(int first, int last)
{
	int i, j, somewhere;
	struct point pivot, temp;



	if (first < last)
	{
		somewhere = (first + last) / 2;
		pivot.x = arr[somewhere].x;
		i = first; j = last;
		do {
			while (arr[i].x < pivot.x) i++;
			while (arr[j].x > pivot.x) j--;
			if (i <= j)
			{
				temp = arr[i]; arr[i] = arr[j]; arr[j] = temp;
				i++; j--;
			}
		} while (i <= j);

		quicksort2(first, j);
		quicksort2(i, last);
	}

}
void quicksort3(int first, int last)
{
	int i, j, somewhere;
	struct point pivot, temp;


	if (first < last)
	{
		somewhere = (first + last) / 2;
		pivot.y = arr[somewhere].y;
		i = first; j = last;
		do {
			while (arr[i].y < pivot.y) i++;
			while (arr[j].y > pivot.y) j--;
			if (i <= j)
			{
				temp = arr[i]; arr[i] = arr[j]; arr[j] = temp;
				i++; j--;
			}
		} while (i <= j);

		quicksort3(first, j);
		quicksort3(i, last);
	}


}
void quicksortO(int first, int last)
{
	int i, j, somewhere;
	struct point pivot, temp;



	if (first < last)
	{
		somewhere = (first + last) / 2;
		pivot.y = inst[somewhere].y;
		i = first; j = last;
		do {
			while (inst[i].y < pivot.y) i++;
			while (inst[j].y > pivot.y) j--;
			if (i <= j)
			{
				temp = inst[i]; inst[i] = inst[j]; inst[j] = temp;
				i++; j--;
			}
		} while (i <= j);

		quicksortO(first, j);
		quicksortO(i, last);
	}

}

int GrahamScan(int m)
{
	int i, top, min;
	struct point temp;



	min = 0;
	for (i = 1; i < m; i++)
		if (p[i].y<p[min].y || p[i].y == p[min].y && p[i].x>p[min].x)
			min = i;
	temp = p[0]; p[0] = p[min]; p[min] = temp;
	quicksort(1, m - 1);
	top = 1;
	for (i = 2; i < m; i++)
	{
		while (area(p[top - 1], p[top], p[i]) <= 0)
			top--;
		top++;
		temp = p[top]; p[top] = p[i]; p[i] = temp;
	}
	return(top + 1);
}

int DelTriangulation()
{
	int n = 0, w = 0;
	int i, j, k, m;
	double xn, yn, zn;
	int flag;
	FILE *fp;

	
	if (!(fp = fopen("createvoronoi3.csv", "r")))//打开EXCEL.csv文件读取数据（A列是x坐标 B列是y坐标）
	{
		printf("无法打开数据文件\n");
		getchar();
		exit(0);

	}
	n = 0;
	while (fscanf(fp, "%lf,%lf", &x[n], &y[n]) != EOF) {
		n++;
	}
	fclose(fp);

	for (i = 0; i < n; i++) {
		z[i] = x[i] * x[i] + y[i] * y[i];

	}

	/*printf("请输入%d个点的x,y坐标", n);
	for (i = 0; i < n; i++)
	{
	scanf("%d %d", &x[i], &y[i]);
	z[i] = x[i] * x[i] + y[i] * y[i];
	}
	for (i = 0; i < n - 2; i++)
		for (j = i + 1; j < n; j++)
			for (k = i + 1; k < n; k++)
				if (j != k)
				{
					xn = (y[j] - y[i])*(z[k] - z[i]) - (y[k] - y[i])*(z[j] - z[i]);
					yn = (x[k] - x[i])*(z[j] - z[i]) - (x[j] - x[i])*(z[k] - z[i]);
					zn = (x[j] - x[i])*(y[k] - y[i]) - (x[k] - x[i])*(y[j] - y[i]);

					if (flag = (zn > 0))
						for (m = 0; m < n; m++) {
							flag = flag && ((x[m] - x[i])*xn + (y[m] - y[i])*yn + (z[m] - z[i])*zn <= 0.00000001);//0.00000001

						}

					if (flag)
					{
						out[w].i = i; out[w].j = j; out[w].k = k;
						printf("out[%d].i=%d out[%d].j=%d out[%d].k=%d\n", w, out[w].i, w, out[w].j, w, out[w].k);
						w++;
					}

				}
	return(w); */
}
double radius3p(double x1, double y1, double x2, double y2, double x3, double y3)
{
	double radius, OX1, OY1;
	OY1 = ((x2*x2*x3) + (x3*y2*y2) - (x1*x1*x3) + (x1*x1*x2) - (y1*y1*x3) + (y1*y1*x2) - (x1*x2*x2) - (x1*y2*y2) - (x2*x3*x3) + (x1*x3*x3) - (x2*y3*y3) + (x1*y3*y3)) / (2 * ((x1*y3) + (y1*x2) + (y2*x3) - (x1*y2) - (x2*y3) - (y1*x3)));
	OX1 = ((x2*x2) + (y2*y2) - (x1*x1) - (y1*y1) + (2 * y1*OY1) - (2 * y2*OY1)) / (2 * (x2 - x1));
	radius = sqrt(pow(x1 - OX1, 2) + pow(y1 - OY1, 2));
	return radius;
}
double finddistance(struct point a, struct point b) {
	double distance;
	distance = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
	return distance;
}
int check_allpoint_that_inside_circle(struct point a, struct point b, int n) {
	double r, d, count = 0, lenght;
	int i;
	struct point o, t;
	d = finddistance(a, b);
	r = d / 2;
	o.x = (a.x + b.x) / 2;
	o.y = (a.y + b.y) / 2;
	for (i = 0; i < n; i++) //先确认凸包上最远两点为直径的圆是否包含所有凸包上的点包含的话就是最小圆
	{
		lenght = finddistance(o, p[i]);
		if (lenght > r) 
		{
			count++;
			printf("最远两点为直径画的圆存在未包含的点,这个点为p[%d]=(%lf,%lf).\n", i + 1, p[i].x, p[i].y);
		}
	}
	if (count == 0) 
	return count;
}


main() {
	int i, j, k,g,l, arrnum = 0, n, vpn = 1, w, check;
	double top = 200, buttom = -200, left = -200, right = 200, y, y2;
	double slope, xmid, ymid, per_bisect, b;
	double slope2, xmid2, ymid2, per_bisect2, b2;
	struct point ins, bij1, bij2, bik1, bik2, c;
	double d, ds, cnt;

	/*3点*/
	double distance, max1l = 0, max2l = 0, max3l = 0;
	int inscount = 0, max1p, max2p, max3p, Oans;
	double radius, minradius = 999;
	struct point temp, max1, max2, max3, ans1, ans2, ans3;
	/*2点*/
	double lenght, maxlenght = 0;
	int outp, m1, m2;
	struct point center;
	double finalradius;
	int ten2 = 0;

	FILE *fp;
	FILE *gp;	// for gnuplot


	if (!(fp = fopen("createvoronoi3.csv", "r")))//打开要代入的EXCEL.csv文件读取数据（A列是x坐标 B列是y坐标）
	{
		printf("无法打开数据文件\n");
		getchar();
		exit(0);

	}
	n = 0;

	while (fscanf(fp, "%lf,%lf", &p[n].x, &p[n].y) != EOF) {
		//printf("p[%d]=(%f %f)\n", n + 1, p[n].x, p[n].y);
		n++;
	}
	fclose(fp);


	/*
	n=10;
	///////乱数/////////////////////
	srand(time(0));
	for (i = 0; i < n; i++) {
		p[i].x = rand() % 100 + 1;
		p[i].y = rand() % 100 + 1;
	}

	for(i=0;i<n;i++){
	printf("%d:",i+1);
	scanf("%d %d",&p[i].x,&p[i].y);
	}
	*/

	

	g = GrahamScan(n);
	printf("---输入的数据坐标---\n");
	for (i = 0; i < n; i++) {
		printf("p[%d]=(%f %f)\n", i + 1, p[i].x, p[i].y);
	}

   w = DelTriangulation();

	for (i = 0; i < g - 1; i++) {
		for (j = i + 1; j < g; j++) {

			slope = (p[j].y - p[i].y) / (p[j].x - p[i].x);

			// 计算中间点

			xmid = (p[i].x + p[j].x) / 2;
			ymid = (p[i].y + p[j].y) / 2;

			// 垂线的斜率

			per_bisect = -1 / (slope);

			// b是在Y轴上的截距
			b = ymid - (per_bisect * xmid);
			y = per_bisect * left + b;
			
			if (per_bisect > 0 && per_bisect<999999) {
				if (y <= buttom) {
					bij1.y = buttom;
					bij1.x = (bij1.y - b) / per_bisect;
					if (per_bisect * right + b >= top) {
						bij2.y = top;
						bij2.x = (top - b) / per_bisect;
					}
					else {
						bij2.x = right;
						bij2.y = per_bisect * bij2.x + b;
					}
				}
				else if (y > buttom) {
					bij1.x = left;
					bij1.y = y;
					if (per_bisect * right + b >= top) {
						bij2.y = top;
						bij2.x = (bij2.y - b) / per_bisect;
					}
					else {
						bij2.x = right;
						bij2.y = per_bisect * bij2.x + b;
					}
				}

			}
			else if (per_bisect < 0 && per_bisect>-999999) {
				if (y <= top) {
					bij1.y = y;
					bij1.x = left;
					if (per_bisect * right + b <= buttom) {
						bij2.y = buttom;
						bij2.x = (bij2.y - b) / per_bisect;
					}
					else {
						bij2.x = right;
						bij2.y = per_bisect * bij2.x + b;
					}
				}
				else {
					bij1.y = top;
					bij1.x = (bij1.y - b) / per_bisect;
					if (per_bisect * right + b <= buttom) {
						bij2.y = buttom;
						bij2.x = (bij2.y - b) / per_bisect;
					}
					else {
						bij2.x = right;
						bij2.y = per_bisect * bij2.x + b;
					}
				}
			}
			else if (per_bisect == 0) {
				bij1.y = b;
				bij1.x = left;
				bij2.y = b;
				bij2.x = right;
			}
			else {
				bij1.y = top;
				bij1.x = xmid;
				bij2.y = buttom;
				bij2.x = xmid;
			}

			arrnum = 0;

			for (k = 0; k < g; k++) {

				if (k != i && k != j) {

					slope2 = (p[k].y - p[i].y) / (p[k].x - p[i].x);

					// 计算中间点

					xmid2 = (p[i].x + p[k].x) / 2;
					ymid2 = (p[i].y + p[k].y) / 2;

					// 垂线的斜率

					per_bisect2 = -1 / (slope2);

					// b是在Y轴上的截距
					b2 = ymid2 - (per_bisect2 * xmid2);
					y2 = per_bisect2 * left + b2;
					
				  if (per_bisect2 > 0 && per_bisect2<999999) {
					
					if (y2 <= buttom) {
						bik1.y = buttom;
						bik1.x = (bik1.y - b2) / per_bisect2;
						if (per_bisect2 *  right + b2 >= top) {
							bik2.y = top;
							bik2.x = (bik2.y - b2) / per_bisect2;
						}
						else {
							bik2.x = right;
							bik2.y = per_bisect2 * bik2.x + b2;
						}
					}
					else if (y2 >= buttom) {
						bik1.x = left;
						bik1.y = y2;
						if (per_bisect2 *  right + b2 >= top) {
							bik2.y = top;
							bik2.x = (bik2.y - b2) / per_bisect2;
						}
						else {
							bik2.x = right;
							bik2.y = per_bisect2 * bik2.x + b2;
						}
					}
					
				   }


				    
					else if (per_bisect2 < 0 && per_bisect2>-999999) {
						if (y2 <= top) {
							bik1.y = y2;
							bik1.x = left;
							if (per_bisect2 * right + b2 <= buttom) {
								bik2.y = buttom;
								bik2.x = (bik2.y - b2) / per_bisect2;
							}
							else {
								bik2.x = right;
								bik2.y = per_bisect2 * bik2.x + b2;
							}
						}
						else {
							bik1.y = top;
							bik1.x = (bik1.y - b2) / per_bisect2;
							if (per_bisect2 * right + b2 <= buttom) {
								bik2.y = buttom;
								bik2.x = (bik2.y - b2) / per_bisect2;
							}
							else {
								bik2.x = right;
								bik2.y = per_bisect2 * bik2.x + b2;
							}
						}
					
					}
					
					
					else if (per_bisect2 == 0) {
						bik1.y = b2;
						bik1.x = left;
						bik2.y = b2;
						bik2.x = right;
					}
					else {
					  bik1.y = top;
					  bik1.x = xmid2;
					  bik2.y = buttom;
					  bik2.x = xmid2;
					}
					ins = intersection(bij1, bij2, bik1, bik2);
					if (ins.x == -999 && ins.y == -999) {
						printf("i=%d j=%d k=%d 这些线段之间没有交叉\n", i, j, k);
					}
					arr[arrnum] = ins;
					arrnum++;
				}
			}
			//put first and last point of bisector line between i and j in intersect array (arr[arrnum]

			if (bij1.x - bij2.x == 0) {
				arr[arrnum] = bij1;
				arrnum++;
				arr[arrnum] = bij2;
				arrnum++;
			}
			else {
				bij1.x = left;
				bij1.y = y;
				bij2.x = right;
				bij2.y = per_bisect * bij2.x + b;
				arr[arrnum] = bij1;
				arrnum++;
				arr[arrnum] = bij2;
				arrnum++;
			}
			check = 1;
			for (l = 1; l <= arrnum - 1; l++) {
				if (arr[l].x == arr[l - 1].x) check++;
			}

			//sort intersect array by x 
			if (check == arrnum) {
				quicksort3(0, arrnum - 1);
			}
			else quicksort2(0, arrnum - 1);

			for (k = 1; k < arrnum; k++) {
				c.x = (arr[k - 1].x + arr[k].x) / 2;
				c.y = (arr[k - 1].y + arr[k].y) / 2;
				d = sqrt(pow(p[i].x - c.x, 2) + pow(p[i].y - c.y, 2));
				cnt = 0;
				for (l = 0; l < g; l++) {

					if (l != i && l != j) {
						ds = sqrt(pow(p[l].x - c.x, 2) + pow(p[l].y - c.y, 2));
						if (ds < d) cnt = cnt + 1;
					}

				}
				if (cnt == g - 2) {
					//printf("(%f,%f)\n(%f,%f)\n", arr[k].x, arr[k].y, arr[k - 1].x, arr[k - 1].y);
					vp[vpn] = arr[k];
					vp[vpn - 1] = arr[k - 1];
					vpn = vpn + 2;
				}
			}
		}
	}
	
	
	



	// gnuplot的启动命令
	if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {	// 通过pipe启动gnuplot
		fprintf(stderr, "没有找到文件 %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	// --- 向gnuplot发送指令 --- //
	fprintf(gp, RANGE);
	fprintf(gp, "set size square\n");
	fprintf(gp, "unset border\n");
	fprintf(gp, "unset xtics\n");
	fprintf(gp, "unset ytics\n");
	fprintf(gp, "unset key\n");
	fprintf(gp, "set style arrow 5 linetype 2 linecolor rgbcolor \'black\' linewidth 4 \n");


	for (i = 1; i < vpn - 1; i = i + 2) {
		fprintf(gp, "set arrow from %lf,%lf to %lf,%lf arrowstyle 5 nohead\n", vp[i].x, vp[i].y, vp[i - 1].x, vp[i - 1].y);
	}


	fprintf(gp, "set term pngcairo enhanced size 1920,1200 transparent truecolor\n");

	fprintf(gp, "set out \"fvoronoi.png\"\n");

	fprintf(gp, "plot'-' with points ls 7 pointsize 1 linecolor rgbcolor \'black\'\n");
	for (i = 0; i < n; i++) {
		fprintf(gp, "%f, %f\n", p[i].x, p[i].y);

	}
	fprintf(gp, "e\n");
	fflush(gp); // 清除缓冲区数据（必须）
	//system("pause");
	fprintf(gp, "exit\n"); // 终止gnuplot
	_pclose(gp);

	
	//////////////////////////////↓↓↓↓↓↓↓↓1中心问题↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓////////////////////////////////
	//////////////////////////////*delete same number in intersection array*///////////////////////////////
	///*
	inscount = 0;
	for (i = 0; i < vpn - 1; i++) {
		if (vp[i].y != buttom && vp[i].y != top && vp[i].x != left && vp[i].x != right) {
			temp = vp[i];
			k = 0;
			for (j = 0; j < inscount; j++) {

				if ((int)temp.x == (int)inst[j].x && (int)temp.y == (int)inst[j].y) {
					k = 1;
				}
			}
			if (k == 0) {

				inst[inscount] = temp;
				inscount++;
			}

		}
	}
	quicksortO(0, inscount - 1);
	printf("---Voronoi图的顶点坐标---\n");
	for (i = 0; i < inscount; i++) {
		printf("O[%d]:（%lf,%lf）\n",i+1, inst[i].x, inst[i].y);
	}

	///////////////////////////////////////////1中心/////////////////////////////////
	/////////////////////////////////////过2点//////////////////////////////////////
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			lenght = finddistance(p[i], p[j]);
			if (lenght >= maxlenght) {
				maxlenght = lenght;
				max1.x = p[i].x; max1.y = p[i].y;
				max2.x = p[j].x; max2.y = p[j].y;
				m1 = i;
				m2 = j;
			}
		}

	}


	outp = check_allpoint_that_inside_circle(max1, max2, n);//确认所有点是不是在圆内



	if (outp == 0) {
		printf("---最小覆盖圆结果---\n");
		printf("点集合中相互距离最远两点的距离为直径的圆包含点集合中所有的点\n");
		printf("两最远点为p[%d]=(%lf,%lf) 和 ", m1 + 1, max1.x, max1.y);//两最远点坐标
		printf("p[%d]=(%lf,%lf)\n", m2 + 1, max2.x, max2.y);
		printf("最小圆直径（最远两点的距离）=%lf\n", maxlenght);
		center.x = (max1.x + max2.x) / 2;
		center.y = (max1.y + max2.y) / 2;
		finalradius = maxlenght / 2;
		ten2 = 1;
	}

	///////////////////////////////////////过3点//////////////////////////////////
	else {

		printf("---Voronoi图各顶点与其最远3母点的坐标以及其最小覆盖圆的半径---\n");
		for (i = 0; i < inscount; i++) {
			max1l = 0; max2l = 0; max3l = 0;

			for (j = 0; j < n; j++) {
				distance = finddistance(inst[i], p[j]);
				if (distance >= max1l && distance >= max2l && distance >= max3l) {
					max1l = distance;
					max1 = p[j];
					max1p = j;
				}
			}
			for (j = 0; j < n; j++) {
				distance = finddistance(inst[i], p[j]);
				if (j != max1p && distance >= max2l && distance >= max3l) {
					max2l = distance;
					max2 = p[j];
					max2p = j;
				}
			}
			for (j = 0; j < n; j++) {
				distance = finddistance(inst[i], p[j]);
				if (j != max1p && j != max2p && distance >= max3l) {
					max3l = distance;
					max3 = p[j];
					max3p = j;
				}
			}
			radius = radius3p(max1.x, max1.y, max2.x, max2.y, max3.x, max3.y);
			if (radius <= minradius) {
				minradius = radius;
				ans1 = max1;
				ans2 = max2;
				ans3 = max3;
				Oans = i;

			}
			printf("O[%d]\n", i + 1);
			printf("坐标:%lf,%lf\n", inst[i].x, inst[i].y);
			printf("点1:%lf,%lf\n", max1.x, max1.y);
			printf("点2:%lf,%lf\n", max2.x, max2.y);
			printf("点3:%lf,%lf\n", max3.x, max3.y);
			radius = radius3p(max1.x, max1.y, max2.x, max2.y, max3.x, max3.y);
			printf("radius=%lf\n\n", radius);
		}
		printf("---最小覆盖圆结果---\n");
		printf("在所有的顶点中，O[%d]形成的圆是最小圆，圆通过的三点是\n", Oans + 1);//最小圆的顶点
		printf("点1:%lf,%lf\n", ans1.x, ans1.y);//最小圆周上的三点坐标（母点）
		printf("点2:%lf,%lf\n", ans2.x, ans2.y);//最小圆周上的三点坐标（母点）
		printf("点3:%lf,%lf\n", ans3.x, ans3.y);//最小圆周上的三点坐标（母点）
		printf("圆半径r=%lf.\n\n", minradius);//最小圆的半径
		center.x = inst[Oans].x;
		center.y = inst[Oans].y;
		finalradius = minradius;

	}
	// gnuplot的启动命令
	if ((gp = _popen(GNUPLOT_PATH, "w")) == NULL) {	// 通过pipe启动gnuplot
		fprintf(stderr, "没有找到文件 %s.", GNUPLOT_PATH);
		exit(EXIT_FAILURE);
	}

	// --- 向gnuplot发送指令 --- //
	fprintf(gp, RANGE);
	fprintf(gp, "set size square\n");
	fprintf(gp, "unset border\n");
	fprintf(gp, "unset xtics\n");
	fprintf(gp, "unset ytics\n");
	fprintf(gp, "unset key\n");
	fprintf(gp, "set style arrow 5 linetype 2 linecolor rgbcolor \'black\' linewidth 4 \n");
	fprintf(gp, "set style arrow 6 linetype 2 linecolor rgbcolor \'black\' linewidth 2 \n");

	for (i = 1; i < vpn - 1; i = i + 2) {
		fprintf(gp, "set arrow from %lf,%lf to %lf,%lf arrowstyle 5 nohead\n", vp[i].x, vp[i].y, vp[i - 1].x, vp[i - 1].y);
	}


	fprintf(gp, "set term pngcairo enhanced size 1920,1200 transparent truecolor\n");

	

	if (ten2 == 1) {
		fprintf(gp, "set arrow from %lf,%lf to %lf,%lf arrowstyle 6 nohead\n", max1.x, max1.y, max2.x, max2.y);
		fprintf(gp, "set out \"fvoronoi1center.png\"\n");
		fprintf(gp, "plot'-' with points ls 7 pointsize 1 linecolor rgbcolor \'black\','-' with points ls 7 pointsize 1 linecolor rgbcolor \'black\', '-' with circle lw 3 linecolor rgbcolor \'black\'\n");
		for (i = 0; i < n; i++) {
			fprintf(gp, "%f, %f\n", p[i].x, p[i].y);

		}
		fprintf(gp, "e\n");
		fprintf(gp, "%f, %f\n", center.x, center.y);
		fprintf(gp, "e\n");
		fprintf(gp, "%lf %lf %lf\n", center.x, center.y, finalradius);
		fprintf(gp, "e\n");
	}
	else {
		fprintf(gp, "set out \"fvoronoi1center.png\"\n");
		fprintf(gp, "plot'-' with points ls 7 pointsize 1 linecolor rgbcolor \'black\', '-' with circle lw 3 linecolor rgbcolor \'black\'\n");
		for (i = 0; i < n; i++) {
			fprintf(gp, "%f, %f\n", p[i].x, p[i].y);

		}
		fprintf(gp, "e\n");
		fprintf(gp, "%lf %lf %lf\n", center.x, center.y, finalradius);
		fprintf(gp, "e\n");

	}

	
	fflush(gp); // 清除缓冲区数据（必须）
	//system("pause");
	fprintf(gp, "exit\n"); // 终止gnuplot
	_pclose(gp);
	//*/
}
