CC=gcc
CFLAGS=-std=c11 -Wall -Wextra
LDFLAGS=-lm
GL_LDFLAGS=-lGL -lglfw

# Files
TARGETS=tiny_ising demo

# Rules
all: $(TARGETS)

tiny_ising: tiny_ising.o ising.o wtime.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

demo: demo.o ising.o wtime.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(GL_LDFLAGS)

clean:
	rm -f $(TARGETS) *.o

.PHONY: clean all



# #File Edit Options Buffers Tools Makefile Help                                                                                                                                                              

# #TODO acordarse de que si modiificamos el makefile y no borramos con clean, no va a hacer nada                                                                                                             
# #TODO si queres automatizar podes hacer esto:                                                                                                                                                              

# #TODO for ((i=0; i<4;i++)); do make clean && make EXTRA_CFLAGS="-O${i}" && ./tiny_mc; done                                                                                                                 
# #TODO otra forma es crear un archivo de parametros.txt                                                                                                                                                     
# #TODO y despues escribir  en la terminal cat parametros.txt | while read p; do make clean && make EXTRA_CFLAGS="${p}" && ./tiny_mc; done                                                                   
# #TODO  while read dicgit re-e cada vez que leas una linea mete lo que dice la linea en la variable p                                                                                                       
# #TODO  se puede cambiar el prompt de la siguiente forma:                                                                                                                                                   
# #TODO cada vez q tiene que imprimir un prompt se fija en una variable de entorno, PS1                                                                                                                      
# #TODO entonces la redefiinis para que valga PS1="#"                                                                                                                                                        
# #TODO  si uno quiere cambiar algo que esdel preprocesador lo pasa con CPPFLAGS, algo asi:                                                                                                                  
# #TODO  do make clean && make CPPFLAGS="-DL=${l}" ahi esta diciendo -DEFINE L                                                                                                                               
# #TODO  https://www.youtube.com/watch?v=YTPXFZIkRV4&t=5113s buen video                                                                                                                                      
# #TODO  se puede cambiar el prompt de la siguiente forma:                                                                                                                                                   

# # Compilers                                                                                                                                                                                                
# CC = gcc

# # Flags                                                                                                                                                                                                    
# EXTRA_CFLAGS=
# CFLAGS = -std=c11 -Wall -Wextra $(EXTRA_CFLAGS)
# LDFLAGS = -lm

# # Binary file                                                                                                                                                                                              
# TARGET = tiny_mc

# # Files                                                                                                                                                                                                    
# C_SOURCES = tiny_mc.c wtime.c
# C_OBJS = $(patsubst %.c, %.o, $(C_SOURCES))

# # Rules se van leyendo las reglas desde arriba hacia abajo                                                                                                                                                 
# all: $(TARGET)

# $(TARGET): $(C_OBJS) #si existe un archivo .o fijate si ves uno .c y si existe compilalo con los parametros de aca abajo                                                                                   
#         $(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
# # $(cc)= gcc                                                                                                                                                                                               
# # $(CFLAGS)= -std=c11 -Wall -Wextra                                                                                                                                                                        
# # $@= es el nombre objetivo de la receta =tiny_mc.o wtime.o                                                                                                                                                
# # $^                                                                                                                                                                                                       
# #       $(CC) $(CFLAGS) -o $@ -c $^ $(LDFLAGS)                                                                                                                                                             

# clean:
#         rm -f $(TARGET) *.o

# .PHONY: clean all



