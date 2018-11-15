#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>

//
// Constantes para OpenGL
//
#define KEY_ESC 27
#define ANCHO 1920
#define ALTO 1080
#define false 0
#define true 1

//
// Constantes para Algoritmo de gravitacion
//
#define PI (3.141592653589793)
#define G 6.673e-11


// ===============
// ===== CPU =====
// ===============

//
// Estructuras y variables para Algoritmo de gravitacion
//

const int CUDA_BLK = 64; // Tamaño predeterminado de bloque de hilos CUDA

float toroide_alfa;
float toroide_theta;
float toroide_incremento;
float toroide_lado;
float toroide_r;
float toroide_R;

int delta_tiempo = 1.0f; //Intervalo de tiempo, longitud de un paso
int pasos;
int N;

//variables nuestras CPU//
float * masas;
float * cPositionX;
float * cPositionY;
float * cPositionZ;

float * cVelocityX;
float * cVelocityY;
float * cVelocityZ;

float * fuerza_totalX;
float * fuerza_totalY;
float * fuerza_totalZ;

//variables nuestras GPU//
float * gpu_masas;
float * gpu_cPositionX;
float * gpu_cPositionY;
float * gpu_cPositionZ;

float * gpu_cVelocityX;
float * gpu_cVelocityY;
float * gpu_cVelocityZ;

float * gpu_fuerza_totalX;
float * gpu_fuerza_totalY;
float * gpu_fuerza_totalZ;

double cColorR = (double) rand() / (RAND_MAX + 1.0);
double cColorG = (double) rand() / (RAND_MAX + 1.0);
double cColorB = (double) rand() / (RAND_MAX + 1.0);

//terminan las variables nuestras.//

//
// Funciones para Algoritmo de gravitacion
//

void calcularFuerzas(int N, int dt) {
    int cuerpo1, cuerpo2;
    float dif_X, dif_Y, dif_Z;
    float distancia;
    float F;

    for (cuerpo1 = 0; cuerpo1 < N - 1; cuerpo1++) {
        for (cuerpo2 = cuerpo1 + 1; cuerpo2 < N; cuerpo2++) {
            if ((cPositionX[cuerpo1] == cPositionX[cuerpo2]) && (cPositionY[cuerpo1] == cPositionY[cuerpo2]) && (cPositionZ[cuerpo1] == cPositionZ[cuerpo2]))
                continue;

            dif_X = cPositionX[cuerpo2] - cPositionX[cuerpo1];
            dif_Y = cPositionY[cuerpo2] - cPositionY[cuerpo1];
            dif_Z = cPositionZ[cuerpo2] - cPositionZ[cuerpo1];

            distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            F = (G * masas[cuerpo1] * masas[cuerpo2]) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            fuerza_totalX[cuerpo1] += dif_X;
            fuerza_totalY[cuerpo1] += dif_Y;
            fuerza_totalZ[cuerpo1] += dif_Z;

            fuerza_totalX[cuerpo2] -= dif_X;
            fuerza_totalY[cuerpo2] -= dif_Y;
            fuerza_totalZ[cuerpo2] -= dif_Z;
        }
    }
}

void moverCuerpos(int N, int dt) {
    int cuerpo;
    for (cuerpo = 0; cuerpo < N; cuerpo++) {

        fuerza_totalX[cuerpo] *= 1 / masas[cuerpo];
        fuerza_totalY[cuerpo] *= 1 / masas[cuerpo];

        cVelocityX[cuerpo] += fuerza_totalX[cuerpo] * dt;
        cVelocityY[cuerpo] += fuerza_totalY[cuerpo] * dt;

        cPositionX[cuerpo] += cVelocityX[cuerpo] * dt;
        cPositionY[cuerpo] += cVelocityY[cuerpo] * dt;

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;

    }
}

void gravitacionCPU(int N, int dt) {
    calcularFuerzas(N, dt);
    moverCuerpos(N, dt);
}

void inicializarEstrella(int i, double n) {

    //todos van a tener la misma masa//
    masas[i] = 0.001 * 8;

    if ((toroide_alfa + toroide_incremento) >= 2 * M_PI) {
        toroide_alfa = 0;
        toroide_theta += toroide_incremento;
    } else {
        toroide_alfa += toroide_incremento;
    }

    cPositionX[i] = (toroide_R + toroide_r * cos(toroide_alfa)) * cos(toroide_theta);
    cPositionY[i] = (toroide_R + toroide_r * cos(toroide_alfa)) * sin(toroide_theta);
    cPositionZ[i] = toroide_r * sin(toroide_alfa);

    cVelocityX[i] = 0.0;
    cVelocityY[i] = 0.0;
    cVelocityZ[i] = 0.0;
}

void inicializarCuerpos(int N) {
    int cuerpo;
    double n = N;

    toroide_alfa = 0.0;
    toroide_theta = 0.0;
    toroide_lado = sqrt(N);
    toroide_incremento = 2 * M_PI / toroide_lado;
    toroide_r = 1.0;
    toroide_R = 2 * toroide_r;

    srand(time(NULL));

    for (cuerpo = 0; cuerpo < N; cuerpo++) {

        fuerza_totalX[cuerpo] = 0.0;
        fuerza_totalY[cuerpo] = 0.0;
        fuerza_totalZ[cuerpo] = 0.0;

        inicializarEstrella(cuerpo, n);

    }
    masas[0] = 2.0e2;
    cPositionX[0] = 0.0;
    cPositionY[0] = 0.0;
    cPositionZ[0] = 0.0;
    cVelocityX[0] = -0.000001;
    cVelocityY[0] = -0.000001;
    cVelocityZ[0] = 0.0;

    masas[1] = 1.0e1;
    cPositionX[1] = -1.0;
    cPositionY[1] = 0.0;
    cPositionZ[1] = 0.0;
    cVelocityX[1] = 0.0;
    cVelocityY[1] = 0.0001;
    cVelocityZ[1] = 0.0;

}

void finalizar(void) {
    free(masas);
    free(cPositionX);
    free(cPositionY);
    free(cPositionZ);
    free(cVelocityX);
    free(cVelocityY);
    free(cVelocityZ);
    free(fuerza_totalX);
    free(fuerza_totalY);
    free(fuerza_totalZ);
}

// ===============
// ===== GPU =====
// ===============

__global__ void kernelCalcularFuerzas(int N, int dt, float * gpu_masas, float * gpu_cPositionX, float * gpu_cPositionY, float * gpu_cPositionZ, float * gpu_cVelocityX, float * gpu_cVelocityY, float * gpu_cVelocityZ, float * gpu_fuerza_totalX, float * gpu_fuerza_totalY, float * gpu_fuerza_totalZ) {

    extern __shared__ float shared_size[];

    //averiguar bien como saltear a la siguiente posición
    float * sh_fuerza_totalX = & shared_size[0];
    float * sh_fuerza_totalY = & shared_size[1 * CUDA_BLK];
    float * sh_fuerza_totalZ = & shared_size[2 * CUDA_BLK];

    float * sh_cPositionX = & shared_size[3 * CUDA_BLK];
    float * sh_cPositionY = & shared_size[4 * CUDA_BLK];
    float * sh_cPositionZ = & shared_size[5 * CUDA_BLK];

    float * sh_cVelocityX = & shared_size[6 * CUDA_BLK];
    float * sh_cVelocityY = & shared_size[7 * CUDA_BLK];
    float * sh_cVelocityZ = & shared_size[8 * CUDA_BLK];

    float * sh_masas = & shared_size[9 * CUDA_BLK];

    //CON ESTA SOLUCIÓN, CONVIENE TENER MÁS HILOS QUE BLOQUES
    //LA MASA SE PUEDE PONER EN MEMORIA CONSTANTE Para eso habría que tenerla en tiempo de ejecución

    //DEFINO EL INDICE DE MI HILO
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int i;
    //Acceden coalescentemente a la memoria global para traer los datos a la memoria compartida
    for (i = 0; i < N; i++) {
        //OBTENGO LOS DATOS DE MASA Y POSICIÓN DE MI INDICE, INICIALIZO VECTOR DE FUERZA EN CERO

        sh_masas[threadIdx.x] = gpu_masas[idx];

        sh_cPositionX[threadIdx.x] = gpu_cPositionX[idx];
        sh_cPositionY[threadIdx.x] = gpu_cPositionY[idx];
        sh_cPositionZ[threadIdx.x] = gpu_cPositionZ[idx];

        sh_cVelocityX[threadIdx.x] = gpu_cVelocityX[idx];
        sh_cVelocityY[threadIdx.x] = gpu_cVelocityY[idx];
        sh_cVelocityZ[threadIdx.x] = gpu_cVelocityZ[idx];
    }

    __syncthreads();
    //Se sincronizan para asegurar que todos los hilos trajeron los datos a procesar

    for (i = 0; i < N; i++) {
        //Fuerzas inicializadas en cero
        sh_fuerza_totalX[i] = 0.0f;
        sh_fuerza_totalY[i] = 0.0f;
        sh_fuerza_totalZ[i] = 0.0f;

        for (int j = 0; j < N; j++) {
            float dif_X = sh_cPositionX[j] - sh_cPositionX[i];
            float dif_Y = sh_cPositionY[j] - sh_cPositionY[i];
            float dif_Z = sh_cPositionZ[j] - sh_cPositionZ[i];

            float distancia = sqrt(dif_X * dif_X + dif_Y * dif_Y + dif_Z * dif_Z);

            float F = (G * sh_masas[i] * sh_masas[j]) / (distancia * distancia);

            dif_X *= F;
            dif_Y *= F;
            dif_Z *= F;

            sh_fuerza_totalX[i] += dif_X;
            sh_fuerza_totalY[i] += dif_Y;
            sh_fuerza_totalZ[i] += dif_Z;
        }
    }

    __syncthreads();
    //Se sincronizan para asegurar que todos los hilos terminaron de procesar

    //Acceden coalescentemente a la memoria global para escribir los resultados
    for (i = 0; i < N; i++) {
        //GUARDO TODOS LOS DATOS EN LA MEMORIA COMPARTIDA, POSICIÓN, FUERZA, VELOCIDAD.
        gpu_fuerza_totalX[idx] = sh_fuerza_totalX[threadIdx.x];
        gpu_fuerza_totalY[idx] = sh_fuerza_totalY[threadIdx.x];
        gpu_fuerza_totalZ[idx] = sh_fuerza_totalZ[threadIdx.x];

        gpu_cPositionX[idx] = sh_cPositionX[threadIdx.x];
        gpu_cPositionY[idx] = sh_cPositionY[threadIdx.x];
        gpu_cPositionZ[idx] = sh_cPositionZ[threadIdx.x];

        gpu_cVelocityX[idx] = sh_cVelocityX[threadIdx.x];
        gpu_cVelocityY[idx] = sh_cVelocityY[threadIdx.x];
        gpu_cVelocityZ[idx] = sh_cVelocityZ[threadIdx.x];
    }
}

__global__ void kernelMoverCuerpos(int N, int dt, float * gpu_masas, float * gpu_cPositionX, float * gpu_cPositionY, float * gpu_cPositionZ, float * gpu_cVelocityX, float * gpu_cVelocityY, float * gpu_cVelocityZ, float * gpu_fuerza_totalX, float * gpu_fuerza_totalY, float * gpu_fuerza_totalZ) {

    //DEFINO EL INDICE DE MI HILO
    /*
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int cuerpo;
    for(cuerpo = 0; cuerpo < N; cuerpo++){

    	gpu_fuerza_totalX[cuerpo] *= 1/gpu_masas[cuerpo];
    	gpu_fuerza_totalY[cuerpo] *= 1/gpu_masas[cuerpo];


    	gpu_cVelocityX[i] += dt * gpu_fuerza_totalX[i];
        gpu_cVelocityY[i] += dt * gpu_fuerza_totalY[i];
        gpu_cVelocityZ[i] += dt * gpu_fuerza_totalZ[i];


    	gpu_cPositionX[cuerpo] += gpu_cVelocityX[cuerpo] *dt;
    	gpu_cPositionY[cuerpo] += gpu_cVelocityY[cuerpo] *dt;


    	gpu_fuerza_totalX[cuerpo] = 0.0;
    	gpu_fuerza_totalY[cuerpo] = 0.0;
    	gpu_fuerza_totalZ[cuerpo] = 0.0;

    }
    */
}

void gravitacionGPU(int N, int dt) {
    //TENEMOS QUE REESCRIBIR ESTAS DOS FUNCIONES PARA QUE ANDEN CON EL GPU

    // Bloque unidimensional de hilos (*blk_size* hilos)
    dim3 dimBlock(CUDA_BLK);

    // Grid unidimensional (*ceil(n/blk_size)* bloques)
    dim3 dimGrid((N + dimBlock.x - 1) / dimBlock.x);

    kernelCalcularFuerzas << < dimGrid, dimBlock >>> (N, dt, gpu_masas, gpu_cPositionX, gpu_cPositionY, gpu_cPositionZ, gpu_cVelocityX, gpu_cVelocityY, gpu_cVelocityZ, gpu_fuerza_totalX, gpu_fuerza_totalY, gpu_fuerza_totalZ);
    cudaDeviceSynchronize();

    //Al kernel de mover cuerpos, le paso todos los parametros de arreglos de la GPU?
    kernelMoverCuerpos << < dimGrid, dimBlock >>> (N, dt, gpu_masas, gpu_cPositionX, gpu_cPositionY, gpu_cPositionZ, gpu_cVelocityX, gpu_cVelocityY, gpu_cVelocityZ, gpu_fuerza_totalX, gpu_fuerza_totalY, gpu_fuerza_totalZ);
    cudaDeviceSynchronize();
}

// ==================
// ===== OpenGL =====
// ==================

//
// Variables OpenGL
//
double alfa = 0.0;

// Para angulo de rotacion y direccion de la camara
float angle = 0.0;
float camAngleX = 0;
float camAngleY = 0;
float distancia = 10;
int ejes = 1;

// Vector actual que representa la direccion de la camara
float lx = 0.0f, lz = -1.0f;
// posicion XZ de la camara
float x = 0.0f, z = 5.0f;

int oldX = 0, oldY = 0;
int rotate = false;

//
// Funciones OpenGL
//

//Funcion que se llama cada vez que se quiere dibujar nuevamente en la pantalla
//Se llama cada vez que se produce el evento render
void GL_camara() {
    float camX, camY, camZ;

    //Camara mirando al origen (pickObjX,pickObjY,pickObjZ) = (0,0,0)
    float pickObjX = 0.0;
    float pickObjY = 0.0;
    float pickObjZ = 0.0;

    camX = distancia * sin(camAngleX);
    camY = distancia * sin(camAngleY);
    camZ = distancia * cos(camAngleY) * cos(camAngleX);

    //Ubicar la camara
    gluLookAt(camX, camY, camZ, // Posicion de la camara
              pickObjX, pickObjY, pickObjZ, // Mirando al punto
              0.0, 1.0, 0.0); // Up vector
}

void GL_dibujarCuerpos(void) {
    int i;

    for (i = 0; i < N; i++) {
        glPushMatrix();
        glTranslatef(cPositionX[i], cPositionY[i], cPositionZ[i]);

        //reemplazar por los valores random de los colores//

        glColor3f((double) rand() / (RAND_MAX + 1.0), (double) rand() / (RAND_MAX + 1.0), (double) rand() / (RAND_MAX + 1.0));
        glutSolidSphere(0.02, 20, 20);

        glPopMatrix();
    }

    //ACA!!! se Llama a la funcion que calcula las fuerzas nuevamente
    //gravitacion GPU//
    gravitacionGPU(N, delta_tiempo);
    //TRAERME LOS DATOS DE LA GPU//

    cudaMemcpy(cPositionX, gpu_cPositionX, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cPositionY, gpu_cPositionY, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cPositionZ, gpu_cPositionZ, N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaMemcpy(cVelocityX, gpu_cVelocityX, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cVelocityY, gpu_cVelocityY, N * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cVelocityZ, gpu_cVelocityZ, N * sizeof(float), cudaMemcpyDeviceToHost);

    //gravitacionCPU(N,delta_tiempo);
}

void GL_dibujar(void) {
    // Borra el color y los buffers de profundidad
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Reiniciar la matriz de transformaciones
    glLoadIdentity();

    //ubica la camara
    GL_camara();

    //Dibuja los ejes de coordenadas (si estan habilitados)
    if (ejes) {
        glBegin(GL_LINES);
        glColor3f(1.0, 0.0, 0.0);
        glVertex3d(0, 0, 0);
        glVertex3d(5.0, 0.0, 0.0);

        glColor3f(0.0, 1.0, 0.0);
        glVertex3d(0, 0, 0);
        glVertex3d(0.0, 5.0, 0.0);

        glColor3f(0.0, 0.0, 1.0);
        glVertex3d(0, 0, 0);
        glVertex3d(0, 0, 5.0);

        glEnd();
    }

    // Dibuja
    glPushMatrix();
    GL_dibujarCuerpos();
    glPopMatrix();

    glutSwapBuffers();
}

void GL_cambioDeDimensionDeVentana(int w, int h) {

    // Evita que se divida por cero cuando la ventana es muy chica
    if (h == 0) h = 1;
    float ratio = w * 1.0 / h;

    // Usa la matriz de proyecion
    glMatrixMode(GL_PROJECTION);

    // Reset matriz
    glLoadIdentity();

    // Configura el viewport para la ventana completa
    glViewport(0, 0, w, h);

    // Configura la perspectiva correcta
    gluPerspective(45.0f, ratio, 0.1f, 100.0f);

    // Modelview
    glMatrixMode(GL_MODELVIEW);
}

//Funcion de inicializacion
void GL_inicio(void) {
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glOrtho(-10, 10, -10, 10, -10, 10);
}

void GL_teclado(unsigned char key, int x, int y) {
    double denominador = 50.0;
    double grados = PI / denominador;
    switch (key) {
    case 'a':
        if (alfa + grados >= 2 * PI)
            alfa = (alfa + grados) - 2 * PI;
        else
            alfa += grados;
        break;
    case '+':
        distancia--;
        break;
    case '-':
        distancia++;
        break;
    case 'e':
        if (ejes == 1) {
            ejes = 0;
        } else {
            ejes = 1;
        }
        break;
    case KEY_ESC:
        finalizar();
        exit(0); // Sale de la aplicacion si se presiona 'Esc'
    }
    glutPostRedisplay();
}

void GL_teclasEspeciales(int key, int x, int y) {
    double denominador = 50.0;
    double grados = PI / denominador;

    switch (key) {
    case GLUT_KEY_RIGHT:
        if (camAngleX - grados < 0)
            camAngleX = (camAngleX - grados) + 2 * PI;
        else
            camAngleX -= grados;
        break;
    case GLUT_KEY_LEFT:
        if (camAngleX + grados >= 2 * PI)
            camAngleX = (camAngleX + grados) - 2 * PI;
        else
            camAngleX += grados;
        break;
    case GLUT_KEY_UP:
        if (camAngleY - grados <= -PI / 2)
            camAngleY = -PI / 2 + 0.001;
        else
            camAngleY -= grados;
        break;
    case GLUT_KEY_DOWN:
        if (camAngleY + grados >= PI / 2)
            camAngleY = PI / 2 - 0.001;
        else
            camAngleY += grados;
        break;
    }

    glutPostRedisplay();
}

void GL_OnMouseDown(int button, int state, int x, int y) {
    rotate = false;
    if (button == GLUT_LEFT_BUTTON) {
        oldX = x;
        oldY = y;
        rotate = true;
    }

}

void GL_OnMouseMove(int x, int y) {

    if (rotate) {
        camAngleX -= (x - oldX) * 0.01f;
        camAngleY += (y - oldY) * 0.01f;
    }

    oldX = x;
    oldY = y;
    glutPostRedisplay();
}

void procesoOpenGL(int argc, char * argv[]) {
    //Inicializa la libreria glut
    glutInit( & argc, argv);
    //Se va a usar doble buffer, paleta RGB
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    //Define la ventana de visualizacion
    glutInitWindowSize(ANCHO, ALTO);

    //Posicionar la ventana
    glutInitWindowPosition(0, 0);
    //Se crea la ventana cuyo nombre en la barra de titulo es lo que viene en argv[0]
    glutCreateWindow(argv[0]);

    //Funcion personalizada que inicializa parametros
    GL_inicio();

    //Define cual es la funcion de control de renderizado
    // Se llama cada vez que se quiere dibujar nuevamente en la pantalla (cada vez que se produce el evento render)
    //GL DIBUJAR LLAMA A NUESTRO CÓDIGO//
    glutDisplayFunc(GL_dibujar);
    glutReshapeFunc(GL_cambioDeDimensionDeVentana);
    glutIdleFunc(GL_dibujar);

    //Define cuales son las funciones que atenderan los eventos del teclado
    glutKeyboardFunc(GL_teclado);
    glutSpecialFunc(GL_teclasEspeciales);

    //Define cuales son las funciones que atenderan los eventos del mouse
    glutMouseFunc(GL_OnMouseDown);
    glutMotionFunc(GL_OnMouseMove);

    //El programa espera aca
    glutMainLoop();
}

int main(int argc, char * argv[]) {

    N = atoi(argv[1]);
    delta_tiempo = atof(argv[2]);
    pasos = atoi(argv[3]);

    //CPU VARIABLES
    cPositionX = (float * ) malloc(N * sizeof(float));
    cPositionY = (float * ) malloc(N * sizeof(float));
    cPositionZ = (float * ) malloc(N * sizeof(float));

    cVelocityX = (float * ) malloc(N * sizeof(float));
    cVelocityY = (float * ) malloc(N * sizeof(float));
    cVelocityZ = (float * ) malloc(N * sizeof(float));

    masas = (float * ) malloc(N * sizeof(float));

    fuerza_totalX = (float * ) malloc(sizeof(float) * N);
    fuerza_totalY = (float * ) malloc(sizeof(float) * N);
    fuerza_totalZ = (float * ) malloc(sizeof(float) * N);

    //GPU VARIABLES
    gpu_cPositionX = (float * ) malloc(N * sizeof(float));
    gpu_cPositionY = (float * ) malloc(N * sizeof(float));
    gpu_cPositionZ = (float * ) malloc(N * sizeof(float));

    gpu_cVelocityX = (float * ) malloc(N * sizeof(float));
    gpu_cVelocityY = (float * ) malloc(N * sizeof(float));
    gpu_cVelocityZ = (float * ) malloc(N * sizeof(float));

    gpu_masas = (float * ) malloc(N * sizeof(float));

    gpu_fuerza_totalX = (float * ) malloc(sizeof(float) * N);
    gpu_fuerza_totalY = (float * ) malloc(sizeof(float) * N);
    gpu_fuerza_totalZ = (float * ) malloc(sizeof(float) * N);

    cudaMalloc( & gpu_cPositionX, N * sizeof(float));
    cudaMalloc( & gpu_cPositionY, N * sizeof(float));
    cudaMalloc( & gpu_cPositionZ, N * sizeof(float));

    cudaMalloc( & gpu_cVelocityX, N * sizeof(float));
    cudaMalloc( & gpu_cVelocityY, N * sizeof(float));
    cudaMalloc( & gpu_cVelocityZ, N * sizeof(float));

    cudaMalloc( & gpu_masas, N * sizeof(float));

    cudaMalloc( & gpu_fuerza_totalX, N * sizeof(float));
    cudaMalloc( & gpu_fuerza_totalY, N * sizeof(float));
    cudaMalloc( & gpu_fuerza_totalZ, N * sizeof(float));

    inicializarCuerpos(N);

    //aca pasamos los datos a la GPU por primera vez//
    cudaMemcpy(gpu_cPositionX, cPositionX, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_cPositionY, cPositionY, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_cPositionZ, cPositionZ, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_cVelocityX, cVelocityX, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_cVelocityY, cVelocityY, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_cVelocityZ, cVelocityZ, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_masas, masas, N * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_fuerza_totalX, fuerza_totalX, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_fuerza_totalY, fuerza_totalY, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_fuerza_totalZ, fuerza_totalZ, N * sizeof(float), cudaMemcpyHostToDevice);

    //ADENTRO DE ESTO SE VA A LLAMAR AL CALCULO//
    procesoOpenGL(argc, argv);

    return (0);

}