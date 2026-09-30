/* New Theory Principles Land Allocation Modeling, an interactive simulation
 *
 * I, the sole author of this software, hereby declare that
 * 	1.) No part of this software, are generated nor derived from large language models,
 *	    and AI-assisted autocomplete technologies.
 *	2.) This software are to be released with MMC051 Math modeling submission.
 *
 * The MIT License
 *
 * Copyright (c) 2026 MMC051 contributor.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
*/

#pragma GCC optimize("Ofast")
#pragma GCC optimize("unroll-loops")
#pragma GCC target("avx2")
#include <bits/stdc++.h>
#include <cassert>
#include "Random.h"
#include "ANSI.h"
using namespace std;
 
#define Point complex<int>
#define PointD complex<float>
#define X real()
#define Y imag()
int Rwater, Rrice, Rcrop;

class Polygon {
	public:
		int W, H;
		int verticesCount;
		vector<Point> vertices;

		void definePolygon(vector<Point> vertices){
			verticesCount = vertices.size();
			this->vertices = vertices;
		}

		int orientationTest(Point P, Point P1, Point P2)
		{
			return (P1.X-P.X)*(P2.Y-P.Y)-(P2.X-P.X)*(P1.Y-P.Y);	
		}
		
		int isInsidePolygon(Point P)
		{
			int count = 0;
			for (int i = 0; i < verticesCount; ++i)
			{
				int nextI = (i + 1) % verticesCount;
				if( (P.Y - vertices[i].Y) * (P.X - vertices[nextI].X) == (P.Y - vertices[nextI].Y) * (P.X - vertices[i].X))
				{
					if( min(vertices[i].X, vertices[nextI].X) <= P.X && P.X <= max(vertices[i].X, vertices[nextI].X)
					&& min(vertices[i].Y, vertices[nextI].Y) <= P.Y && P.Y <= max(vertices[i].Y, vertices[nextI].Y) )
					{
						// Point is located on the boundary
						return 1;
						break;
					}
				}
				if (i == verticesCount - 1)
				{
					int count = 0;
					for (int j = 0; j < verticesCount; ++j)
					{
						int nextJ = (j + 1) % verticesCount;
						if(vertices[j].Y <= P.Y && P.Y < vertices[nextJ].Y && orientationTest(P, vertices[j], vertices[nextJ])>0) ++count;
						if(vertices[nextJ].Y <= P.Y && P.Y < vertices[j].Y && orientationTest(P, vertices[nextJ], vertices[j])>0) ++count;
					}
					return count % 2;	
				}
			}
			return 0;
		}
};


class Utilities
{
	public:
		inline float toDegrees(float x){
			return x*numbers::pi/180;
		}

		inline float dist(int ch1,int ch2,vector<PointD> vec){
			return sqrt((vec[ch1].X-vec[ch2].X)*(vec[ch1].X-vec[ch2].X)+(vec[ch1].Y-vec[ch2].Y)*(vec[ch1].Y-vec[ch2].Y));
		}
};

class ZonePolygon: public Polygon, public Utilities {
	public:
		int grid[107][007];
		int zone[107][107];
		
		int Z, theta, gamma;
		float sinTheta, cosTheta, tanGamma;
		float dmax=0, hmax=0;
		vector<Point> latticePoints;
		
		vector<PointD> centerOfMass;
		vector<float> centerOfHeight;
		vector<float> Ai;

		inline float h(PointD P){
			return P.X * sinTheta * tanGamma + P.Y * cosTheta * tanGamma;
		}

		void zoneDivision()
		{
			int count = 0;
			for (int x = 0; x < W; ++x)		
			{
				for (int y = 0; y < H; ++y)
				{
					int location = isInsidePolygon(Point(x, y));
					int rng=Random::get(1,100);
					if(location)
					{
						grid[x][y] = location;
						latticePoints.push_back(Point(x,y));
						hmax = max(hmax, h(PointD(x,y)));
						
						if(rng<=Rwater) rng=1;
						else if(rng<=Rwater+Rrice) rng=2;
						else if(rng<=Rwater+Rrice+Rcrop) rng=3;
						else rng=4;
						
						zone[x][y]=rng;
						Ai[rng]++;
						count++;
					}
				}
			}
			
			for(int i=0;i<latticePoints.size();++i)
			{
				/*
				int rng;
				if(i<0.4*latticePoints.size()) rng = 2;
				else if(i<0.7*latticePoints.size()) rng = 1;
				else if(i<0.9*latticePoints.size()) rng = 3;
				else rng = 4;
				zone[latticePoints[i].X][latticePoints[i].Y]=rng;
				*/
			}

			
			for(auto i:vertices){
				for(auto j:vertices){
					dmax=max(dmax,(float)sqrt((i.X-j.X)*(i.X-j.X)+(i.Y-j.Y)*(i.Y-j.Y)));
				}
			}
		}
		
		int depth[107][107];	
		float V;
		inline void updateCenter(){
			fill(Ai.begin(),Ai.end(),0);
			queue<pair<pair<int,int>,int>> q;
			vector<float> sx,sy,sh;
			sx.resize(Z+1); sy.resize(Z+1); sh.resize(Z+1);
			for(int i=0;i<W;++i){
				for(int j=0;j<H;++j){
					depth[i][j]=0;
					if (zone[i][j] == 1)
					{
						if(i == 0 || i== W-1 || j == 0 || j == H-1){
							depth[i][j] = 1;
						}		
						else {	
							for (int r1 = -1; r1 <= 1; ++r1){
								for (int r2 = -1; r2 <=1 ; ++r2){
									if(zone[i+r1][j+r2]!=zone[i][j])
									{
										depth[i][j]=1;
									}
								}
							}
						}
						if(depth[i][j])
						{
							q.push({{i,j},1});
						}
					}
					sx[zone[i][j]]+=i;
					sy[zone[i][j]]+=j;
					sh[zone[i][j]]+=h(PointD(i,j));
					++Ai[zone[i][j]];		
				}
			}
			for(int i=1;i<=Z;++i){
				sx[i]/=Ai[i]; sy[i]/=Ai[i]; sh[i]/=Ai[i];
				centerOfMass[i]={sx[i],sy[i]};
				centerOfHeight[i]=sh[i];
			
			}
			while(!q.empty()){
				int t1 = q.front().first.first;
				int t2 = q.front().first.second;
				int t3 = q.front().second;
				q.pop();
				for (int r1 = -1; r1 <= 1; ++r1){
					for(int r2 = -1; r2 <= 1; ++r2){
						if(t1+r1>=W || t2+r2>=H || t1+r1<0 || t2+r2<0) continue;
						if(zone[t1+r1][t2+r2]!=1) continue;
						if(!depth[t1 + r1][t2 + r2]){
							depth[t1 + r1][t2 + r2] = t3 + 1;
							q.push({{t1 + r1, t2 + r2}, t3 + 1});
						}
					}
				}
			}
			V = 0;
			for(int i=0;i<W;++i){
				for(int j=0;j<H;++j){
					V += depth[i][j]/3.0;
				}
			}
		}
		
};
class Land: public ZonePolygon
{
	public:
		int DR, DW, Rin;
		float U2 = 0.0090, U3 = 0.0045;
		float dayOpt = 365, kW1 = 7.4511e-6;
		float Wd2 = 0.3119, Wh2 = 0.4905, Wl2 = 0.1976;
		float Wd3 = 0.5485, Wh3 = 0.2409, Wl3 = 0.2106, kh3wohmax = -log(0.5);
		// Initialization
		Land(vector<Point> vertices, int theta, int gamma, int DR, int DW, int Rin){	
			definePolygon(vertices);
			this->W = 50;
			this->H = 80;
			this->Z = 4;
			this->theta = theta;
			this->gamma = gamma;
			this->DR = DR;
			this->DW = DW;
			this->Rin = Rin;
			sinTheta = sin(toDegrees(theta));
			cosTheta = cos(toDegrees(theta));
			tanGamma = tan(toDegrees(gamma));
			centerOfMass.resize(Z+1);
			centerOfHeight.resize(Z+1);
			Ai.resize(Z+1);
			zoneDivision();
			updateCenter();
		}
	
		// Print

		inline void print(){
			ANSI::reset();
			for(int i=0;i<W;++i){
				for(int j=0;j<H;++j){
					if(zone[i][j]==1) ANSI::color(34);
					if(zone[i][j]==2) ANSI::color(33);
					if(zone[i][j]==3) ANSI::color(32);
					if(zone[i][j]==4) ANSI::color(31);
					//ANSI::color(30+zone[i][j]);
					if(zone[i][j]==0) cout << ' ';
					else cout << zone[i][j];
				}
				cout << endl;
			}
			ANSI::color(0);
			updateCenter();	
			cout << "STATISTIC " << totalEfficiency(centerOfMass, centerOfHeight, V, Ai) << " "  << endl;
			for(auto i:Ai) cout << i << " ";
			cout << "W VALIDITY " << (V>=60*U2*Ai[2]+60*U3*Ai[3]);
			/*for(int i=0;i<W;++i){
				for(int j=0;j<H;++j) cout << depth[i][j];
				cout << endl;
			}*/
		}

		// Optimizer
		float prev = 0;
		int counter = 0;
		inline void optimizer(int time){
			updateCenter();
			auto oldc1 = centerOfMass;
			auto oldc2 = centerOfHeight; 
			auto oldc3 = V;
			auto oldc4 = Ai;
			

			vector<int>  a(time),b(time),ax(time),ay(time),bx(time),by(time),var(time);
			for(int i=0;i<time;++i){
				a[i] = Random::get(0,latticePoints.size()-1);
				b[i] = Random::get(0,latticePoints.size()-1);
				ax[i]=latticePoints[a[i]].X; ay[i]=latticePoints[a[i]].Y;
				bx[i]=latticePoints[b[i]].X; by[i]=latticePoints[b[i]].Y;
				
				if(zone[ax[i]][ay[i]] !=3 ){
					do
					{
						int dx = (ax[i]-centerOfMass[zone[ax[i]][ay[i]]].X>0)? -1:1;
						int dy = (ay[i]-centerOfMass[zone[ax[i]][ay[i]]].Y>0)? -1:1;
						dx *= Random::get(0,20);
						dy *= Random::get(0,20);
						bx[i]=ax[i]+dx,by[i]=ay[i]+dy;
					} while(bx[i]<0 || by[i]<0 || bx[i]>=W || by[i]>=H || find(latticePoints.begin(),latticePoints.end(),Point(bx[i],by[i]))==latticePoints.end());
				}
				var[i]=Random::get(1,Z);
				swap(zone[ax[i]][ay[i]],zone[bx[i]][by[i]]);
			}
			
			
			///swap(zone[ax][ay],var);
			updateCenter();
			auto newc1 = centerOfMass;
			auto newc2 = centerOfHeight;
			auto newc3 = V;
			auto newc4 = Ai;
			// Constraints
			//
			
			if(totalEfficiency(newc1, newc2, newc3, newc4)>=totalEfficiency(oldc1, oldc2, oldc3, oldc4)){
				// Favorable
				prev = V;
				counter = 0;
			}
			else{
				for(int i=time-1;i>=0;--i){
					swap(zone[ax[i]][ay[i]],zone[bx[i]][by[i]]);
				}
			
				/*
				//swap(zone[ax][ay],zone[bx][by]);
				if(counter>10){
					swap(zone[ax][ay],zone[bx][by]);
					//swap(zone[ax][ay],var);
				}
				else{
					prev=V;
				}
				++counter;
				*/
			}

		}

		
		// Model

		inline float totalEfficiency(vector<PointD> centerOfMass, vector<float> centerOfHeight, float V, vector<float> Ai){
			if(!(Rin>=U2*Ai[2]+U3*Ai[3])) return -1; 				// 3.5.1 Water mechanics
			//if(!(V>=60*U2*Ai[2]+60*U3*Ai[3]) && V<prev) return -1;		// 3.5.2 Volume	
			float Atotal = Ai[1]+Ai[2]+Ai[3]+Ai[4];
			if(!(Ai[2]/Atotal>=300.00/1600.00)) return -1;
			if(!(Ai[3]>=200.00/1600.00)) return -1;
			if(!(Ai[4]>=100.00/1600.00)) return -1;
			float statistic = 0;
			for (int zone = 1; zone <= Z; ++zone)
			{
				float Wi, Ei;
				
				if(zone == 1)
				{
					float dayActual = V / (U2 * Ai[2] + U3 * Ai[3]);
					Ei = exp(-kW1*(dayActual - dayOpt)*(dayActual - dayOpt));
					Wi = 0.3416;
				}
				else if(zone == 2){
					float d21=dist(2,1,centerOfMass);
					float Ed2 = 1 - abs(d21/dmax);
					float Eh2 = 1 - centerOfHeight[2]/hmax;
					float phi2 = 90 + atan(toDegrees((centerOfMass[4].Y-centerOfMass[2].Y)/(centerOfMass[4].X-centerOfMass[2].X)));
					float El2 = phi2/180;
					Ei = (Ed2*Wd2+Eh2*Wh2+El2*Wl2)/(Wd2+Wh2+Wl2);
					Wi = 0.3774;
				}
				else if(zone == 3){
					float d31=dist(3,1,centerOfMass);
					float Ed3 = 1 - abs(d31/dmax);
					float Eh3 = exp(-kh3wohmax/(hmax*hmax/4.0)*(centerOfHeight[3]-(hmax/2))*(centerOfHeight[3]-(hmax/2)));
					float phi3 = 90 + atan(toDegrees((centerOfMass[4].Y-centerOfMass[3].Y)/(centerOfMass[4].X-centerOfMass[3].X)));
					float El3 = phi3/180;
					Ei = (Ed3*Wd3+Eh3*Wh3+El3*Wl3)/(Wd3+Wh3+Wl3);
					//cout << d31 << " " << Ed3 << " " << Eh3 << " " << phi3 << " " << El3 << " " << Ei << endl;	
					Wi = 0.1405;
				}
				else if(zone == 4){
					Ei = centerOfHeight[4]/hmax;
					Wi = 0.1258;
				}
				else{
					Ei = 0;
				}
				//cout << Ei << " ";
				statistic += Ai[zone] * Wi * Ei;
			}
			//cout << statistic  << " " << endl;
			return statistic;
		}


};


signed main(){
	cout << "Input theta, gamma, iterations: ";
	int theta, gamma, iterations, temp1, temp2; cin >> theta >> gamma >> iterations;
	cout << "Input Rwater, Rrice, Rcrop, Rresidential (sum = 100)";
	cin >> Rwater >> Rrice >> Rcrop >> temp1;
	int N; cin >> N;
	vector<Point> input;
	for (int i = 0; i < N; ++i){
		cin >> temp1 >> temp2;
		input.push_back(Point(temp1, temp2));
	}	
	Land farm(input,theta,gamma,0,0,300);
	farm.print();
	for(int i=0;i<iterations;++i){
		farm.optimizer(2);
		if(i%1000==0) farm.print();
	}
	farm.print();
		
}
