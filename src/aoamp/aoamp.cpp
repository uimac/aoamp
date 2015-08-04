/*

aobench C code is licensed under 2-clause BSD.

Copyright 2009-2014, Syoyo Fujita
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#include <amp.h>
#include <amp_math.h>
using namespace concurrency;

#include "amp_tinymt_rng.h"

#include <mmsystem.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#define M_PI 3.1415
#define WIDTH        256
#define HEIGHT       256
#define NSUBSAMPLES  2
#define NAO_SAMPLES  8

float drand48()
{
	return ((float)(rand()) / RAND_MAX);
}

typedef struct _vec
{
	float x;
	float y;
	float z;
} vec;


typedef struct _Isect
{
	float t;
	vec    p;
	vec    n;
	int    hit;
} Isect;

typedef struct _Sphere
{
	vec    center;
	float radius;

} Sphere;

typedef struct _Plane
{
	vec    p;
	vec    n;

} Plane;

typedef struct _Ray
{
	vec    org;
	vec    dir;
} Ray;

Sphere spheres[3];
Plane  plane;

static float vdot(vec v0, vec v1) restrict(amp)
{
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

static void vcross(vec *c, vec v0, vec v1) restrict(amp)
{

	c->x = v0.y * v1.z - v0.z * v1.y;
	c->y = v0.z * v1.x - v0.x * v1.z;
	c->z = v0.x * v1.y - v0.y * v1.x;
}

static void vnormalize(vec *c) restrict(amp)
{
	float length = concurrency::fast_math::sqrt(vdot((*c), (*c)));

	if (concurrency::fast_math::fabs(length) > 1.0e-17) {
		c->x /= length;
		c->y /= length;
		c->z /= length;
	}
}

void
ray_sphere_intersect(Isect *isect, const Ray *ray, const Sphere *sphere) restrict(amp)
{
	vec rs;

	rs.x = ray->org.x - sphere->center.x;
	rs.y = ray->org.y - sphere->center.y;
	rs.z = ray->org.z - sphere->center.z;

	float B = vdot(rs, ray->dir);
	float C = vdot(rs, rs) - sphere->radius * sphere->radius;
	float D = B * B - C;

	if (D > 0.0) {
		float t = -B - concurrency::fast_math::sqrt(D);

		if ((t > 0.0) && (t < isect->t)) {
			isect->t = t;
			isect->hit = 1;

			isect->p.x = ray->org.x + ray->dir.x * t;
			isect->p.y = ray->org.y + ray->dir.y * t;
			isect->p.z = ray->org.z + ray->dir.z * t;

			isect->n.x = isect->p.x - sphere->center.x;
			isect->n.y = isect->p.y - sphere->center.y;
			isect->n.z = isect->p.z - sphere->center.z;

			vnormalize(&(isect->n));
		}
	}
}

void
ray_plane_intersect(Isect *isect, const Ray *ray, const Plane *plane) restrict(amp)
{
	float d = -vdot(plane->p, plane->n);
	float v = vdot(ray->dir, plane->n);

	if (concurrency::fast_math::fabs(v) < 1.0e-4) return;

	float t = -(vdot(ray->org, plane->n) + d) / v;

	if ((t > 0.0) && (t < isect->t)) {
		isect->t = t;
		isect->hit = 1;

		isect->p.x = ray->org.x + ray->dir.x * t;
		isect->p.y = ray->org.y + ray->dir.y * t;
		isect->p.z = ray->org.z + ray->dir.z * t;

		isect->n = plane->n;
	}
}

void
orthoBasis(vec *basis, vec n) restrict(amp)
{
	basis[2] = n;
	basis[1].x = 0.0; basis[1].y = 0.0; basis[1].z = 0.0;

	if ((n.x < 0.6) && (n.x > -0.6)) {
		basis[1].x = 1.0;
	}
	else if ((n.y < 0.6) && (n.y > -0.6)) {
		basis[1].y = 1.0;
	}
	else if ((n.z < 0.6) && (n.z > -0.6)) {
		basis[1].z = 1.0;
	}
	else {
		basis[1].x = 1.0;
	}

	vcross(&basis[0], basis[1], basis[2]);
	vnormalize(&basis[0]);

	vcross(&basis[1], basis[2], basis[0]);
	vnormalize(&basis[1]);
}


void ambient_occlusion(vec *col, const Isect *isect, tinymt& rand, array_view<Sphere, 1> sphere_view, array_view<Plane, 1> plane_view) restrict(amp)
{
	int    i, j;
	int    ntheta = NAO_SAMPLES;
	int    nphi = NAO_SAMPLES;
	float eps = 0.0001f;

	vec p;

	p.x = isect->p.x + eps * isect->n.x;
	p.y = isect->p.y + eps * isect->n.y;
	p.z = isect->p.z + eps * isect->n.z;

	vec basis[3];
	orthoBasis(basis, isect->n);

	float occlusion = 0.0;

	for (j = 0; j < ntheta; j++) {
		for (i = 0; i < nphi; i++) {
			float theta = concurrency::fast_math::sqrt(rand.next_single());
			float phi = static_cast<float>(2.0f * M_PI * rand.next_single());

			float x = concurrency::fast_math::cos(phi) * theta;
			float y = concurrency::fast_math::sin(phi) * theta;
			float z = concurrency::fast_math::sqrt(1.0f - theta * theta);

			// local -> global
			float rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
			float ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
			float rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;

			Ray ray;

			ray.org = p;
			ray.dir.x = rx;
			ray.dir.y = ry;
			ray.dir.z = rz;

			Isect occIsect;
			occIsect.t = static_cast<float>(1.0e+17);
			occIsect.hit = 0;

			ray_sphere_intersect(&occIsect, &ray, &sphere_view[0]);
			ray_sphere_intersect(&occIsect, &ray, &sphere_view[1]);
			ray_sphere_intersect(&occIsect, &ray, &sphere_view[2]);
			ray_plane_intersect(&occIsect, &ray, &plane_view[0]);

			if (occIsect.hit) occlusion += 1.0;

		}
	}

	occlusion = (ntheta * nphi - occlusion) / (float)(ntheta * nphi);

	col->x = occlusion;
	col->y = occlusion;
	col->z = occlusion;
}

unsigned char
clamp(float f)
{
	int i = (int)(f * 255.5);

	if (i < 0) i = 0;
	if (i > 255) i = 255;

	return (unsigned char)i;
}


void
render(unsigned char *img, int w, int h, int nsubsamples)
{
	float *fimg = (float *)malloc(sizeof(float) * w * h * 3);
	memset((void *)fimg, 0, sizeof(float) * w * h * 3);

	// fimgをGPUに送る.
	array_view<float, 1> fimg_view(w * h * 3, fimg);
	fimg_view.discard_data();

	// sphereをGPUに送る.
	array_view<Sphere, 1> sphere_view(3, spheres);

	// planeをGPUに送る
	array_view<Plane, 1> plane_view(1, &plane);

	// indexをGPUに送る.
	std::vector<int> index_data(w * h);
	for (int i = 0, wh = w * h; i < wh; ++i) {
		index_data[i] = i;
	}
	array_view<int, 1> index_view(w * h, index_data);

	// 乱数.
	tinymt_collection<1> myrand(extent<1>(w * h), rand());

	parallel_for_each(
		index_view.extent,
		[=](index<1> idx) restrict(amp)
	{
		const int pos = index_view[idx];
		const int x = pos % (h);
		const int y = pos / (w);
		const index<1> p0(pos * 3 + 0);
		const index<1> p1(pos * 3 + 1);
		const index<1> p2(pos * 3 + 2);
		tinymt& rand = myrand[idx];

		for (int v = 0; v < nsubsamples; v++) {
			for (int u = 0; u < nsubsamples; u++) {
				float px = static_cast<float>((x + (u / (float)nsubsamples) - (w / 2.0)) / (w / 2.0));
				float py = static_cast<float>(-(y + (v / (float)nsubsamples) - (h / 2.0)) / (h / 2.0));

				Ray ray;

				ray.org.x = 0.0;
				ray.org.y = 0.0;
				ray.org.z = 0.0;

				ray.dir.x = px;
				ray.dir.y = py;
				ray.dir.z = -1.0;
				vnormalize(&(ray.dir));

				Isect isect;
				isect.t = static_cast<float>(1.0e+17);
				isect.hit = 0;

				ray_sphere_intersect(&isect, &ray, &sphere_view[0]);
				ray_sphere_intersect(&isect, &ray, &sphere_view[1]);
				ray_sphere_intersect(&isect, &ray, &sphere_view[2]);
				ray_plane_intersect(&isect, &ray, &plane_view[0]);

				if (isect.hit) {
					vec col;
					ambient_occlusion(&col, &isect, rand, sphere_view, plane_view);

					fimg_view[p0] += col.x;
					fimg_view[p1] += col.y;
					fimg_view[p2] += col.z;
				}

			}
		}

		fimg_view[p0] /= (float)(nsubsamples * nsubsamples);
		fimg_view[p1] /= (float)(nsubsamples * nsubsamples);
		fimg_view[p2] /= (float)(nsubsamples * nsubsamples);
	}
	);

	for (int i = 0, size = w * h * 3; i < size; ++i) {
		img[i] = clamp(fimg_view[i]);
	}
}

void
init_scene()
{
	spheres[0].center.x = -2.0f;
	spheres[0].center.y = 0.0f;
	spheres[0].center.z = -3.5f;
	spheres[0].radius = 0.5f;

	spheres[1].center.x = -0.5f;
	spheres[1].center.y = 0.0f;
	spheres[1].center.z = -3.0f;
	spheres[1].radius = 0.5f;

	spheres[2].center.x = 1.0f;
	spheres[2].center.y = 0.0f;
	spheres[2].center.z = -2.2f;
	spheres[2].radius = 0.5f;

	plane.p.x = 0.0f;
	plane.p.y = -0.5f;
	plane.p.z = 0.0f;

	plane.n.x = 0.0f;
	plane.n.y = 1.0f;
	plane.n.z = 0.0f;

}

void
saveppm(const char *fname, int w, int h, unsigned char *img)
{
	FILE *fp;

	fp = fopen(fname, "wb");
	assert(fp);

	fprintf(fp, "P6\n");
	fprintf(fp, "%d %d\n", w, h);
	fprintf(fp, "255\n");
	fwrite(img, w * h * 3, 1, fp);
	fclose(fp);
}

int
main(int argc, char **argv)
{
	DWORD time = ::timeGetTime();
	unsigned char *img = (unsigned char *)malloc(WIDTH * HEIGHT * 3);

	init_scene();

	render(img, WIDTH, HEIGHT, NSUBSAMPLES);

	printf("%f \n", (::timeGetTime() - time) / 1000.0);

	saveppm("ao.ppm", WIDTH, HEIGHT, img);

	return 0;
}
