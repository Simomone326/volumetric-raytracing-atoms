#include <SDL3/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#define WIDTH 800
#define HEIGHT 600
#define LENGTH 100
#define STEPNUM 200
#define COLORPAR 15000
#define SCALE 0.05

#define A0 1
#ifndef PI
#define PI 3.14159
#endif

//this saves the image
void saveBufferToBMP(uint32_t* buffer, int width, int height, const char* filename) {
    // 1. Create a surface that points to your existing pixel buffer
    // Pitch is the width of one row in bytes (width * 4 bytes for uint32)
    SDL_Surface* surface = SDL_CreateSurfaceFrom(
        width, 
        height, 
        SDL_PIXELFORMAT_ARGB8888, 
        buffer, 
        width * sizeof(uint32_t)
    );

    if (surface == NULL) {
        printf("Failed to create surface for saving: %s\n", SDL_GetError());
        return;
    }

    // 2. Save the surface to a BMP file
    if (!SDL_SaveBMP(surface, filename)) {
        printf("Failed to save BMP: %s\n", SDL_GetError());
    } else {
        printf("Successfully saved orbital to %s\n", filename);
    }

    // 3. Clean up the surface (this does NOT free your buffer, just the surface wrapper)
    SDL_DestroySurface(surface);
}

typedef struct{
    float x;
    float y;
    float z;
} Vector3;


double assoc_legendre(int l, int m, double x) {
    // 1. Domain checks
    if (m < 0 || m > l || fabs(x) > 1.0) return 0.0;

    // 2. Compute P_m^m(x)
    // Formula: (-1)^m * (2m-1)!! * (1-x^2)^(m/2)
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    
    if (l == m) return pmm;

    // 3. Compute P_{m+1}^m(x)
    // Formula: x * (2m + 1) * P_m^m
    double pmmp1 = x * (2.0 * m + 1.0) * pmm;
    
    if (l == m + 1) return pmmp1;

    // 4. Compute P_l^m(x) using recurrence
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll = (x * (2.0 * ll - 1.0) * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
        
        pmm = pmmp1;
        pmmp1 = pll;
    }

    return pll;
}

double assoc_laguerre(int n, int alpha, double x) {
    if (n < 0) return 0.0;
    if (n == 0) return 1.0;

    // L_0^alpha(x) = 1
    double L_prev = 1.0; 
    
    // L_1^alpha(x) = 1 + alpha - x
    double L_curr = 1.0 + alpha - x;

    if (n == 1) return L_curr;

    // Recurrence for n >= 2
    for (int k = 2; k <= n; k++) {
        double L_next = ((2 * k - 1 + alpha - x) * L_curr - (k - 1 + alpha) * L_prev) / k;
        
        L_prev = L_curr;
        L_curr = L_next;
    }

    return L_curr;
}

double factorial(int n){
    double ret = 1;
    for(int i = n; i > 0; i--){
        ret*=i;
    }
    return ret;
}

float sqrtInRadialComp(int l,int n){
    float ret = (2.0/(n*A0))*(2.0/(n*A0))*(2.0/(n*A0)) * (((double)factorial(n-l-1))/(2.0*n*(double)factorial(n+l)));
    ret = sqrtf(ret);
    return ret;
}

float sqrtInAngularComp(int l, int m){
    float ret = ((2.0*l+1)/(4*PI))*(((double)factorial(l-m))/((double)factorial(l+m)));
    ret = sqrtf(ret);
    return ret;
}

float radialComp(int n, int l, float r, float radical){
    return radical*exp(-(r/(n*A0)))*pow((2.0*r)/(n*A0), l)*assoc_laguerre(n-l-1, 2*l+1, (2.0*r)/(n*A0));   
}

float angularComp(float radical,int l,int m,float theta,float phi){
    float P_val = assoc_legendre(l, abs(m), cos(theta));

    // 2. Add the Phase / Directionality (The "Real" mapping)
    if (m == 0) {
        return radical*P_val;
    } 
    else if (m > 0) {
        // "Real" component (e.g., aligns with X-axis)
        // Usually requires a normalization factor of sqrt(2), but for shape it doesn't matter
        return radical * sqrt(2.0) * P_val * cos(m * phi);
    } 
    else { // m < 0
        // "Imaginary" component (e.g., aligns with Y-axis)
        return radical * sqrt(2.0) * P_val * sin(abs(m) * phi);
    }
}


uint32_t makeColor(int r, int g, int b) {
    // 1. Clamp values to make sure they stay between 0 and 255
    if (r > 255) r = 255; if (r < 0) r = 0;
    if (g > 255) g = 255; if (g < 0) g = 0;
    if (b > 255) b = 255; if (b < 0) b = 0;

    uint32_t a = 255; // We want the pixel to be fully visible (opaque)

    // 2. Shift the bits into their correct "slots"
    // Alpha moves 24 bits left, Red moves 16, Green moves 8, Blue stays put.
    return (a << 24) | (r << 16) | (g << 8) | b;
}

Vector3 newpos(Vector3 prev, Vector3 u, float step){
    Vector3 ret;
    ret.x = prev.x + u.x*step;
    ret.y = prev.y + u.y*step;
    ret.z = prev.z + u.z*step;
    return ret;
}


Vector3 vecSub(Vector3 v1, Vector3 v2){
    Vector3 ret = {v1.x - v2.x, v1.y - v2.y, v1.z - v2.z};
    return ret;
}

float mod(Vector3 v){
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

void normalize(Vector3 *pV){
    float _mod = mod(*pV);
    (*pV).x/=_mod;
    (*pV).y/=_mod;
    (*pV).z/=_mod;
}


float calcWaveFunc(int n,int l,int m, float r, float theta, float phi, float angRadical,float radRadical){
    return radialComp(n, l, r, radRadical)*angularComp(angRadical, l, m, theta, phi);
}

void calcPos(float *r, float *phi, float *theta, Vector3 pos, Vector3 center){

    Vector3 sub_v = {pos.x - center.x, pos.y - center.y, pos.z - center.z};

    sub_v.x *= SCALE;
    sub_v.y *= SCALE;
    sub_v.z *= SCALE;

    float r_val = mod(sub_v);
    //division by 0
    if (r_val < 1e-6) r_val = 1e-10;
    if (fabs(sub_v.x) < 1e-6) sub_v.x = 1e-10;

    float cosTheta = sub_v.z/r_val;
    //for float imprecision
    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta < -1.0) cosTheta = -1.0;

    

    *r = r_val;
    *phi = atan2(sub_v.y, sub_v.x);
    *theta = acos(cosTheta);
    


}

float calculateV(Vector3 pos, Vector3 atomPos, int N, int L, int M){

    float r, theta, phi;
    calcPos(&r, &phi, &theta, pos, atomPos);
    
    float angRadical = sqrtInAngularComp(L, M);
    float radialRadical = sqrtInRadialComp(L, N);

    float PSI = calcWaveFunc(N, L, M, r, theta, phi, angRadical, radialRadical);
    float PSIsq = PSI*PSI;
    return PSIsq;
}

void getPixels(uint32_t buffer[], Vector3 position, int N, int L, int M){
    

    Vector3 camera = {position.x, position.y-50, position.z};
    const Vector3 CENTER = {0, 10, 0};
    Vector3 startPxlPos = {position.x - WIDTH/2, position.y, position.z + HEIGHT/2};
    for (int i = 0; i < WIDTH*HEIGHT; i++)
    {
        //idxPxl = (y*HEIGHT) + x
        int pxlX = i % WIDTH;
        int pxlY = (i - pxlX) / WIDTH;
        Vector3 pxlPos = {startPxlPos.x + pxlX, startPxlPos.y, startPxlPos.z - pxlY};
        Vector3 u = vecSub(pxlPos, camera);
        normalize(&u);
        Vector3 ray = pxlPos;
        float v = 0;
        for(int j = 0; j < STEPNUM; j++){
            ray = newpos(ray, u, (float)LENGTH/STEPNUM);

            
            v += calculateV(ray, CENTER, N, L, M)*COLORPAR;
        }      

        buffer[i] = makeColor(v, 0, v);
    }

}




int main(int argc, char *argv[]) {

    //read n, l, m
    int N, L, M;
    char save;
    printf("write n, l and m separated by a space:\n");
    scanf("%d %d %d", &N, &L, &M);
    printf("n:%d l:%d m:%d\n",N, L, M);
    printf("do you want to save the image? (Y/N): ");
    scanf(" %c", &save);

    // 1. Setup SDL
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("Error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;
    SDL_CreateWindowAndRenderer("Atom", WIDTH, HEIGHT, 0, &window, &renderer);

    // 2. Create the "Canvas" (Texture)
    // ARGB8888 means: Alpha (Transparency), Red, Green, Blue (8 bits each)
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);

    // 3. Create the pixel buffer in standard C memory
    // We use uint32_t so 1 number = 1 pixel (4 bytes)
    uint32_t *buffer = malloc(WIDTH * HEIGHT * sizeof(uint32_t));

    int running = 1;
    SDL_Event event;

    
    Vector3 position = {0, 0, 0};
    for (int i = 0; i < WIDTH * HEIGHT; i++) {
            buffer[i] = 0xFF000000; 
        }
    getPixels(buffer, position, N, L, M);
    
    // --- Event Loop ---
    while (running) {

        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) running = 0;
        }
  
        // 4. Update the Texture with our new pixel data
        // "WIDTH * sizeof(uint32_t)" is the "pitch" (bytes per row)
        SDL_UpdateTexture(texture, NULL, buffer, WIDTH * sizeof(uint32_t));

        // 5. Render
        SDL_RenderClear(renderer);
        SDL_RenderTexture(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }
    
    
    //save
    
    if(save == 'Y'){
        char filename[50];
        sprintf(filename, "orbital n=%d l=%d m=%d.bmp", N, L, M);
        saveBufferToBMP(buffer, WIDTH, HEIGHT, filename);
    }

    // Cleanup
    free(buffer);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    printf("over");

    return 0;
}