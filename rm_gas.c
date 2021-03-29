/************************************************
Title: rm_gas.c
Author: Jared Coughlin
Date: 10/27/15
Purpose: 2LPTIC seems incapable of making ics with
         only dm, so this code simply removes the
         gas particles.
Notes:
************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>



/***********************
      Globals
***********************/
int HEADER_SIZE = 256;
char snap_filename[200];
int NumPartTot;
int num_withmassesTot;
typedef struct IO_HEADER
{
	int npart[6];
	double mass[6];
	double time;
	double redshift;
   int flag_sfr;
	int flag_feedback;
	int npartTotal[6];
	int flag_cooling;
	int num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int flag_stellarage;
	int flag_metals;
	unsigned int npartTotalHighWord[6];
	int flag_entropy_instead_u;
   int flag_RII;
   int flag_RIa;
	char flag_Carbon;
	char flag_Nitrogen;
	char flag_Oxygen;
	char flag_Florine;
	char flag_Neon;
	char flag_Sodium;
	char flag_Magnesium;
	char flag_Aluminum;
	char flag_Silicon;
	char flag_Phosphorus;
	char flag_Sulfur;
	char flag_Chlorine;
	char flag_Argon;
	char flag_Potassium;
	char flag_Calcium;
	char flag_Scandium;
	char flag_Titanium;
	char flag_Vanadium;
	char flag_Chromium;
	char flag_Manganese;
	char flag_Iron;
	char flag_Cobalt;
	char flag_Nickel;
	char flag_Copper;
	char flag_Zinc;
	char fill[27];
} IO_HEADER;

typedef struct PARTICLE_DATA
{
	float Pos[3];
	float Vel[3];
	int Type, id;
	float Mass, U;
} PARTICLE_DATA;



// Fields for block checking
enum fields
{
   HEADER,
   POS,
   VEL,
   IDS,
   MASS,
   U,
};


PARTICLE_DATA *M_P;
IO_HEADER m_header;



/***********************
      Prototypes
***********************/
void single_file_read(void);
IO_HEADER get_header(char *);
int get_particles(IO_HEADER);
int get_withmass(IO_HEADER);
PARTICLE_DATA *read_snapshot(char *, int, int, IO_HEADER);
size_t my_fread(void *, size_t, size_t, FILE *);
size_t my_fwrite(void *, size_t, size_t, FILE *);
void block_check(enum fields, int, int, IO_HEADER);
int get_block_size(enum fields, IO_HEADER);
void rm_gas(char *);



/***********************
         main
***********************/
int main(int argc, char **argv)
{
   // Check args
   if(argc != 2)
   {
      printf("Error, incorrect number of args! Just need snapshot name,\n");
      printf("e.g. snapshot_001\n");
      exit(EXIT_FAILURE);
   }

   // Read in the args
   sprintf(snap_filename, "./%s", argv[1]);

   // Read the snapshot
   single_file_read();

   // Remove gas particles
   rm_gas(argv[1]);

   return 0;
}



/***********************
   single_file_read
***********************/
void single_file_read(void)
{
	/* This is the driver routine for reading in the data from 
 		 the snapshot if there is only 1 file per snapshot. */

	// First we read in the header
	m_header = get_header(snap_filename);

	// Now get the number of particles
	NumPartTot = get_particles(m_header);
	
	// Now get the number of particles with different masses
	num_withmassesTot = get_withmass(m_header);
	
	// Now read in the the data
	M_P = read_snapshot(snap_filename, NumPartTot, num_withmassesTot, m_header);
}



/***********************
      get_header
***********************/
IO_HEADER get_header(char *filename)
{
	/* This function reads in the snapshot header and returns it. */

	IO_HEADER header;
	FILE *fd;
	int dummy;
	#define SKIP my_fread(&dummy, sizeof(dummy), 1, fd);

	// Open file for reading.
	if(!(fd = fopen(filename, "rb")))
	{
		printf("Error, could not open file for reading!\n");
		exit(EXIT_FAILURE);
	}

	// Read header
	SKIP;
	my_fread(&header, sizeof(IO_HEADER), 1, fd);
	fclose(fd);

	return header;
}



/***********************
     get_particles
***********************/
int get_particles(IO_HEADER header)
{
	/* This function returns the number of particles in the 
 		 current file being read. */

	int i, NumPart;

	for(i = 0, NumPart = 0; i < 6; i++)
	{
		NumPart += header.npart[i];
	}

	return NumPart;
}



/***********************
     get_withmass
***********************/
int get_withmass(IO_HEADER header)
{
	/* This function returns the number of particles that 
 		 have differing masses (and are thus present in the mass 
		 block). */

	int i, n_withmass;

	for(i = 0, n_withmass = 0; i < 6; i++)
	{
		if(header.mass[i] == 0)
		{
			n_withmass += header.npart[i];
		}
	}

	return n_withmass;
}



/***********************
    read_snapshot
***********************/
PARTICLE_DATA *read_snapshot(char* fullname, int NumPart, int ntot_withmasses, IO_HEADER header1)
{
	/* This function actually reads in the data from the current file. */

	FILE *fd;
	int n, k, blksize1, blksize2, pc = 1, pc_new, pc_sph; //pc_star;
   enum fields blocknr;
	PARTICLE_DATA *P;

	if(!(P = malloc(NumPart * sizeof(PARTICLE_DATA))))
	{
		printf("Could not allocate memory for particles.\n");
		exit(0);
	}	
	P--;

	IO_HEADER dummy_header;
	#define SKIP1 my_fread(&blksize1, sizeof(blksize1), 1, fd);
   #define SKIP2 my_fread(&blksize2, sizeof(blksize2), 1, fd);

	/*Open file*/
	if(!(fd = fopen(fullname, "rb")))
	{
		printf("Can't open '%s' for reading snapshot.\n", fullname);
		exit(0);
	}

	/*Read dummy header*/
	SKIP1;
	my_fread(&dummy_header, sizeof(dummy_header), 1, fd);
	SKIP2;
   
   blocknr = HEADER;
   block_check(blocknr, blksize1, blksize2, dummy_header);

	/*Positions and type*/
	SKIP1;
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < header1.npart[k]; n++)
		{
			my_fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
			P[pc_new].Type = k;								
			pc_new++;
		}
	}
	SKIP2;

   blocknr = POS;
   block_check(blocknr, blksize1, blksize2, dummy_header);

	/*Velocities*/
	SKIP1;
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < header1.npart[k]; n++)
		{
			my_fread(&P[pc_new].Vel[0], sizeof(float), 3, fd);
			pc_new++;
		}
	}
	SKIP2;

   blocknr = VEL;
   block_check(blocknr, blksize1, blksize2, dummy_header);

	/*Ids*/
	SKIP1;
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < header1.npart[k]; n++)
		{
			my_fread(&P[pc_new].id, sizeof(int), 1, fd);														
			pc_new++;
		}
	}
	SKIP2;

   blocknr = IDS;
   block_check(blocknr, blksize1, blksize2, dummy_header);

	/*Masses*/
	if(ntot_withmasses > 0)
	{
		SKIP1;
	}
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < header1.npart[k]; n++)
		{
			P[pc_new].Type = k;
			if(header1.mass[k] == 0)
			{
				my_fread(&P[pc_new].Mass,sizeof(float),1, fd);
			}
			else
			{
				P[pc_new].Mass = header1.mass[k];
			}
			pc_new++;
		}
	}
	if(ntot_withmasses > 0)
	{
		SKIP2;

      blocknr = MASS;
      block_check(blocknr, blksize1, blksize2, dummy_header);
	}
	
	/*Clean up*/
	fclose(fd);

	return P;
}



/***********************
        my_fread
***********************/
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
   /* This function is taken from gadget. It checks to make sure fread has 
      read the correct number of elements. It handles the compiler warning about 
      ignoring the return value of the function when using -O3. */

   size_t nread;

   if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
   {
      printf("I/O error with fread.\n");
      exit(EXIT_FAILURE);
   }

   return nread;
}



/***********************
        my_fwrite
***********************/
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
   /* This function is taken from gadget. It checks to make sure fwrite has 
      written the correct number of elements. It handles the compiler warning about 
      ignoring the return value of the function when using -O3. */

   size_t nwritten;

   if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
   {
      printf("I/O error with fwrite.\n");
      exit(EXIT_FAILURE);
   }

   return nwritten;
}



/***********************
      block_check
***********************/
void block_check(enum fields blocknr, int blksize1, int blksize2, IO_HEADER dummy_header)
{
   // This function gets the size of block and makes sure that
   // the padding values that are written around it (the ones 
   // obtained with SKIP) match that value.

   int size;

   // Make sure the two padding values match
   if(blksize1 != blksize2)
   {
      printf("Paddings don't match! Block: %d\n", blocknr);
      exit(EXIT_FAILURE);
   }

   // Make sure the padding matches the actual size of the block
   size = get_block_size(blocknr, dummy_header);

   if(blksize1 != size)
   {
      printf("Paddings do not match actual block size! Block: %d\n", blocknr);
      exit(EXIT_FAILURE);
   }
}



/***********************
    get_block_size
***********************/
int get_block_size(enum fields blocknr, IO_HEADER h)
{
   // This function gets the actual block size to make sure that the correct
   // values are written to the padding.

   int bsize;

   switch(blocknr)
   {
      case HEADER:
         bsize = sizeof(IO_HEADER);
         break;

      case POS:
      case VEL:
         bsize = sizeof(float) * 3 * (h.npart[0] + h.npart[1] + h.npart[2] + 
                 h.npart[3] + h.npart[4] + h.npart[5]); 
         break;

      case IDS:
         bsize = sizeof(int) * (h.npart[0] + h.npart[1] + h.npart[2] + h.npart[3] +
                 h.npart[4] + h.npart[5]);
         break;

      case MASS:
         bsize = sizeof(float) * get_withmass(h);
         break;

      case U:
         bsize = sizeof(float) * h.npart[0];
         break;
   }

   return bsize;
}



/***********************
         rm_gas
***********************/
void rm_gas(char *filename)
{
   // Removes gas particles from M_P. By this, rather than fucking with
   // memory and such, it just rewrites the file but skips the gas particles
   // and changes the appropriate variables in the header

	FILE *fd;
	int n, k, blksize, pc = 1, pc_new;
   int ntot_withmasses;
   int ngas;
   char newname[200];

	#define SKIP3 my_fwrite(&blksize, sizeof(blksize), 1, fd);

	/*Open file*/
   sprintf(newname, "%s.dm", filename);

	if(!(fd = fopen(newname, "wb")))
	{
		printf("Can't open '%s' for writing snapshot.\n", newname);
		exit(0);
	}

   // Save number of gas particles for skipping purposes
   ngas = m_header.npart[0];

   // Change the number of gas particles to 0 in header
   m_header.npart[0] = 0;
   m_header.npartTotal[0] = 0;

	/*Write header*/
   blksize = sizeof(m_header);

	SKIP3;
	my_fwrite(&m_header, sizeof(m_header), 1, fd);
	SKIP3;   

	/*Positions*/
   blksize = (m_header.npart[0] + m_header.npart[1] + m_header.npart[2] + 
              m_header.npart[3] + m_header.npart[4] + 
              m_header.npart[5]) * sizeof(float) * 3.0;
	SKIP3;
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < m_header.npart[k]; n++)
		{
			my_fwrite(&M_P[ngas + pc_new].Pos[0], sizeof(float), 3, fd);								
			pc_new++;
		}
	}
	SKIP3;

	/*Velocities*/
	SKIP3;
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < m_header.npart[k]; n++)
		{
			my_fwrite(&M_P[ngas + pc_new].Vel[0], sizeof(float), 3, fd);
			pc_new++;
		}
	}
	SKIP3;

	/*Ids*/
   blksize = (m_header.npart[0] + m_header.npart[1] + m_header.npart[2] + 
              m_header.npart[3] + m_header.npart[4] + 
              m_header.npart[5]) * sizeof(int);
	SKIP3;
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < m_header.npart[k]; n++)
		{
			my_fwrite(&M_P[ngas + pc_new].id, sizeof(int), 1, fd);														
			pc_new++;
		}
	}
	SKIP3;

	/*Masses*/
   ntot_withmasses = get_withmass(m_header);
   blksize = ntot_withmasses * sizeof(float);

	if(ntot_withmasses > 0)
	{
		SKIP3;
	}
	for(k = 0, pc_new = pc; k < 6; k++)
	{
		for(n = 0; n < m_header.npart[k]; n++)
		{
			if(m_header.mass[k] == 0)
			{
				my_fwrite(&M_P[ngas + pc_new].Mass,sizeof(float),1, fd);
			}
			pc_new++;
		}
	}
	if(ntot_withmasses > 0)
	{
		SKIP3;
	}
	
	/*Clean up*/
	fclose(fd);
}
