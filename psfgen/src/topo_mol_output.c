/***************************************************************************
 *cr
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *      $RCSfile: topo_mol_output.c,v $
 *      $Author: jribeiro $        $Locker:  $             $State: Exp $
 *      $Revision: 1.46 $      $Date: 2020/03/10 04:54:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *  
 ***************************************************************************/

#include <string.h>
#include <stdlib.h>
#include "topo_mol_output.h"
#include "topo_mol_struct.h"
#include "pdb_file.h"
#include "psfgen.h"

int topo_mol_write_pdb(topo_mol *mol, FILE *file, void *vdata, void *v, 
                                void (*print_msg)(void *, void *, const char *)) {

  char buf[128], insertion[2];
  int iseg,nseg,ires,nres,atomid,resid;
  int has_guessed_atoms = 0;
  double x,y,z,o,b;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
    
#if defined(NEWPSFGEN)

  int i;

#endif

  
  if ( ! mol ) return -1;

  write_pdb_remark(file,"original generated coordinate pdb file");

  atomid = 0;
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;

    if ( strlen(seg->segid) > 4 ) {
      sprintf(buf,
	"warning: truncating segid %s to 4 characters allowed by pdb format",
	seg->segid);
      print_msg(vdata, v,buf);
    }

    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)  
      for ( atom = res->atoms; atom; atom = atom->next ) {
#else
      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
#endif
        /* Paranoid: make sure x,y,z,o are set. */
        x = y = z = 0.0; o = -1.0;
        ++atomid;
        
        switch ( atom->xyz_state ) {
        case TOPO_MOL_XYZ_SET:
          x = atom->x;  y = atom->y;  z = atom->z;  o = 1.0;
          break;
        case TOPO_MOL_XYZ_GUESS:
        case TOPO_MOL_XYZ_BADGUESS:
          x = atom->x;  y = atom->y;  z = atom->z;  o = 0.0;
          has_guessed_atoms = 1;
          break;
        default:
          print_msg(vdata, v,"ERROR: Internal error, atom has invalid state.");
          print_msg(vdata, v,"ERROR: Treating as void.");
          /* Yes, fall through */
        case TOPO_MOL_XYZ_VOID:
          x = y = z = 0.0;  o = -1.0;
          break;
        }
        b = atom->partition;
        insertion[0] = 0;
        insertion[1] = 0;
        sscanf(res->resid, "%d%c", &resid, insertion);


#if !defined(NOIO)
        write_pdb_atom(file,atomid,atom->name,res->name,resid,insertion,
		(float)x,(float)y,(float)z,(float)o,(float)b,res->chain,
		seg->segid,atom->element);
#endif

#if defined(NEWPSFGEN)
        /* Add the drude particle */
        if (atom->alpha) {
          ++atomid;
#if !defined(NOIO)
          write_pdb_atom(file,atomid,atom->dname,res->name,resid,insertion,
  		atom->dxyz[0],atom->dxyz[1],atom->dxyz[2],(float)o,(float)b,res->chain,
  		seg->segid,atom->element);
#endif      
        }
#endif


      }
    }
  }

  write_pdb_end(file);
  if (has_guessed_atoms) {
    print_msg(vdata, v, 
        "Info: Atoms with guessed coordinates will have occupancy of 0.0.");
  }
  return 0;
}

int topo_mol_write_namdbin(topo_mol *mol, FILE *file, FILE *velfile, void *vdata,
                   void *v, void (*print_msg)(void *, void *, const char *)) {

  int numatoms, iseg, nseg, ires, nres, has_void_atoms=0;
  double x, y, z, xyz[3];
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;

#if defined(NEWPSFGEN)
  int i;
#endif

  if ( ! mol ) return -1;

#if defined(NEWPSFGEN)
  /* stop the writing process if there are drude
   * particles until the namd binary format is updated.
   * Since the lonepairs are treated as atoms at the pdb level, and there is
   * no extra field to add to the file (alpha and thole), we can use
   * this command to write structures with lonepairs
   */
  if (mol->drude) {
    print_msg(vdata, v, "ERROR: NAMD binary file is not yet compatible with Drude" 
      " force-field. Please use the commands writepsf/" 
      "writepdb.");
    return -1;
  }
#endif

  numatoms = 0;
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)  
      for ( atom = res->atoms; atom; atom = atom->next ) {
        ++numatoms;
      }
#else
      for ( i = 0; i < res->atomSize; i++ ) {
        ++numatoms;
      }
#endif

    }
  }
  if ( fwrite(&numatoms, sizeof(int), 1, file) != 1 ) {
    print_msg(vdata, v, "error writing namdbin file");
    return -2;
  }

  if ( velfile ) {
    if ( fwrite(&numatoms, sizeof(int), 1, velfile) != 1 ) {
      print_msg(vdata, v, "error writing velnamdbin file");
      return -4;
    }
  }

  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)  
      for ( atom = res->atoms; atom; atom = atom->next ) {
#else
      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
#endif

        /* Paranoid: make sure x,y,z are set. */
        x = y = z = 0.0;
        switch ( atom->xyz_state ) {
        case TOPO_MOL_XYZ_SET:
        case TOPO_MOL_XYZ_GUESS:
        case TOPO_MOL_XYZ_BADGUESS:
          x = atom->x;  y = atom->y;  z = atom->z;
          break;
        default:
          print_msg(vdata, v,"ERROR: Internal error, atom has invalid state.");
          print_msg(vdata, v,"ERROR: Treating as void.");
          /* Yes, fall through */
        case TOPO_MOL_XYZ_VOID:
          x = y = z = 0.0;
          has_void_atoms = 1;
          break;
        }
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
        if ( fwrite(xyz, sizeof(double), 3, file) != 3 ) {
          print_msg(vdata, v, "error writing namdbin file");
          return -3;
        }

#if defined(NEWPSFGEN)
        /* Add the drude particle */
        if (atom->alpha) {
          if ( fwrite(xyz, sizeof(double), 3, file) != 3 ) {
            print_msg(vdata, v, "error writing namdbin file");
            return -3;
          }
        }
#endif

        if ( velfile ) {
          xyz[0] = atom->vx;
          xyz[1] = atom->vy;
          xyz[2] = atom->vz;
          if ( fwrite(xyz, sizeof(double), 3, velfile) != 3 ) {
            print_msg(vdata, v, "error writing velnamdbin file");
            return -5;
          }

#if defined(NEWPSFGEN)
          /* Add the drude particle */
          if (atom->alpha) {
            if ( fwrite(xyz, sizeof(double), 3, file) != 3 ) {
              print_msg(vdata, v, "error writing namdbin file");
              return -3;
            }
          }
#endif

        }
      }
    }
  }

  if (has_void_atoms) {
    print_msg(vdata, v, "Warning: Atoms with unknown coordinates written at 0. 0. 0.");
  }

  return 0;
}


int topo_mol_write_psf(topo_mol *mol, FILE *file, int charmmfmt, int nocmap, 
                       int nopatches, void *vdata, void *v, 
                       void (*print_msg)(void *, void *, const char *)) {
  char buf[128];
  int iseg,nseg,ires,nres,atomid;
  int namdfmt, charmmext;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  topo_mol_bond_t *bond;
  int nbonds;
  topo_mol_angle_t *angl;
  int nangls;
  topo_mol_dihedral_t *dihe;
  int ndihes;
  topo_mol_improper_t *impr;
  int nimprs;
  topo_mol_cmap_t *cmap;
  int ncmaps;
  topo_mol_exclusion_t *excl;
  int nexcls;
  int numinline;
  int npres,ipres,ntopo,itopo;
  topo_defs_topofile_t *topo;
  topo_mol_patch_t *patch;
  topo_mol_patchres_t *patchres;
  char defpatch[10];
  fpos_t ntitle_pos, save_pos;
  char ntitle_fmt[128];
  int ntitle_count;

#if defined(NEWPSFGEN)
  int i, j;
  double tmpalpha, tmpthole, tmpmass, tmpcharge;
  int numlonepairs, numlphosts, numdrude, numaniso;
  int atomuint; /* unknown int coloumn before the alpha and thole (drude ff),
                 * or last column in the atom section that is -1 in the case
                 * of lonepairs
                 */
  topo_mol_anisotropy_t *aniso;
  float lpunit; /* units to be applied to the angles and dihedrals of the 
                  * bisector and relative lone pair (180.0/M_PI). 
                  * Colinear should be 1.0
                  */
  // Get the psfge_data structure from data
  psfgen_data *psfcontext = *(psfgen_data **)vdata; 
#endif
  strcpy(defpatch,"");

  if ( ! mol ) return -1;

  namdfmt = 0;
  charmmext = 0;
  atomid = 0;
  nbonds = 0;
  nangls = 0;
  ndihes = 0;
  nimprs = 0;
  ncmaps = 0;
  nexcls = 0;
  
#if defined(NEWPSFGEN)
  numdrude = 0;
  numlonepairs = 0;
  numlphosts = 0;
  numaniso = 0;
  lpunit = 1.0;
#endif
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;

    if ( strlen(seg->segid) > 4 ) {
      charmmext = 1;
    }

    nres = hasharray_count(seg->residue_hash);
    
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      if (strlen(res->resid) > 4) {
        charmmext = 1;
      }

#if !defined(NEWPSFGEN)  

      for ( atom = res->atoms; atom; atom = atom->next ) {
        atom->atomid = ++atomid;
#else
      /* count the anisotropy entries */
      if (res->numaniso) {
        numaniso += res->numaniso;
      }
      /* update the atomids and search for drude hosts 
       * atom->alpha != 0. This has to be done prior to the 
       * definition of bonds and etc. and the drude particles
       * increment the atomid with 1. All the atomdis need to be updated 
       * accordingly
      */
      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        atom->atomid = ++atomid;
        if (atom->alpha) {
          ++atomid;
          ++numdrude;
          
          if (psfcontext->VPBONDS) ++nbonds;
          
        }
      }
      
      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
         if (atom->isdrudlonepair) {
          ++numlonepairs;
          switch (atom->lonepair->lptype) {
            case COLINEARLP: 
              numlphosts += 2;
              break;
            case RELATIVELP:
            case BISECTORLP:
              numlphosts += 3;
              break;
            default:
              numlphosts += 3;
              break;
          }
          if (psfcontext->VPBONDS) ++nbonds;
        }
#endif

        

#if !defined(NEWPSFGEN)
        if (strlen(atom->name) > 4) {
          charmmext = 1;
        }
        if ((! charmmfmt) && (strlen(atom->type) > 4)) {
          charmmext = 1;
        }
        if ((! charmmfmt) && (strlen(atom->type) > 6)) {
          namdfmt = 1;
        }
#else
      if (!charmmext) {
        if (strlen(atom->name) > 4) {
          charmmext = 1;
        }
        if ((! charmmfmt) && (strlen(atom->type) > 4)) {
          charmmext = 1;
        }
        
        if ((! namdfmt) && (strlen(atom->type) > 6)) {
          namdfmt = 1;
        }
      }
#endif  

#if !defined(NEWPSFGEN)
        for ( bond = atom->bonds; bond;
                bond = topo_mol_bond_next(bond,atom) ) {
          if ( bond->atom[0] == atom && ! bond->del ) {
            ++nbonds;
          }
        }
            
#else

        
        for ( bond = atom->bonds; bond;
                bond = topo_mol_bond_next(bond,atom) ) {
          if ( bond->atom[0]->atomid == atom->atomid && ! bond->del) {
            
            if (bond->atom[0]->isdrudlonepair || 
                bond->atom[1]->isdrudlonepair ) {
              continue;
            }
            ++nbonds;
          }
        
        }
#endif       
        for ( angl = atom->angles; angl;
                angl = topo_mol_angle_next(angl,atom) ) {
          if ( angl->atom[0] == atom && ! angl->del ) {
            ++nangls;
          }
        }
        
#if  !defined(NEWPSFGEN)        
        for ( dihe = atom->dihedrals; dihe;
                dihe = topo_mol_dihedral_next(dihe,atom) ) {
          if ( dihe->atom[0] == atom && ! dihe->del ) {
            ++ndihes;
          }
        }
#else 
        for ( dihe = atom->dihedrals; dihe;
                dihe = dihe->next ) {
          if (!dihe->del && 
              !dihe->atom[0]->del && 
              !dihe->atom[1]->del && 
              !dihe->atom[2]->del) ++ndihes;
        }
#endif

#if  !defined(NEWPSFGEN) 
        for ( impr = atom->impropers; impr;
                impr = topo_mol_improper_next(impr,atom) ) {
          if ( impr->atom[0] == atom && ! impr->del ) {
            ++nimprs;
          }
        }
#else
      for ( impr = atom->impropers; impr;
              impr = impr->next ) {
        if (!impr->del && 
            !impr->atom[0]->del && 
            !impr->atom[1]->del && 
            !impr->atom[2]->del) {
          ++nimprs;
        }
      }
#endif
        for ( cmap = atom->cmaps; cmap;
                cmap = topo_mol_cmap_next(cmap,atom) ) {
          if ( cmap->atom[0] == atom && ! cmap->del ) {
            ++ncmaps;
          }
        }
        for ( excl = atom->exclusions; excl;
                excl = topo_mol_exclusion_next(excl,atom) ) {
          if ( excl->atom[0] == atom && ! excl->del ) {
            ++nexcls;
          }
        }
      }
    }
  }
  sprintf(buf,"total of %d atoms",atomid);
  print_msg(vdata, v,buf);
  
#if defined(NEWPSFGEN)
  if (numdrude) {
    sprintf(buf,"total of %d drude particles",numdrude);
    print_msg(vdata, v,buf);
  }
  if (numaniso) {
    sprintf(buf,"total of %d anisotropy entries",numaniso);
    print_msg(vdata, v,buf);
  }
  if (numlonepairs) {
    sprintf(buf,"total of %d lone pairs",numlonepairs);
    print_msg(vdata, v,buf);
  }
#endif

  sprintf(buf,"total of %d bonds",nbonds);
  print_msg(vdata, v,buf);
  sprintf(buf,"total of %d angles",nangls);
  print_msg(vdata, v,buf);
  sprintf(buf,"total of %d dihedrals",ndihes);
  print_msg(vdata, v,buf);
  sprintf(buf,"total of %d impropers",nimprs);
  print_msg(vdata, v,buf);
  sprintf(buf,"total of %d explicit exclusions",nexcls);
  print_msg(vdata, v,buf);

  if ( namdfmt ) { charmmext = 0; }
  else if ( atomid > 9999999 ) { charmmext = 1; }

  if ( namdfmt ) {
    print_msg(vdata, v,"Structure requires space-delimited NAMD PSF format");
  } else if ( charmmext ) {
    print_msg(vdata, v,"Structure requires EXTended PSF format");
  }

  ntitle_fmt[0] = '\0';
  strcat(ntitle_fmt, "PSF");
  if ( namdfmt ) strcat(ntitle_fmt, " NAMD");
  if ( charmmext ) strcat(ntitle_fmt, " EXT");
  if ( nocmap ) {
    sprintf(buf,"total of %d cross-terms (not written to file)",ncmaps);
  } else {
    sprintf(buf,"total of %d cross-terms",ncmaps);
    if ( ncmaps ) {
      strcat(ntitle_fmt, " CMAP");
    } else {
      nocmap = 1;
    }
  }
  #if defined(NEWPSFGEN)
    if ( numdrude ) strcat(ntitle_fmt, " DRUDE");
  #endif
  print_msg(vdata, v,buf);
  strcat(ntitle_fmt, charmmext ? "\n\n%10d !NTITLE\n" : "\n\n%8d !NTITLE\n");

#if !defined(NOIO)
  fgetpos(file,&ntitle_pos);
  fprintf(file,ntitle_fmt,1);
#endif
  ntitle_count = 1;

#if !defined(NOIO)
  if ( charmmfmt ) 
   fprintf(file," REMARKS %s\n","original generated structure charmm psf file");
  else
   fprintf(file," REMARKS %s\n","original generated structure x-plor psf file");
#endif
  
  if (mol->npatch) {
    ntitle_count++;
    
#if !defined(NOIO)
    fprintf(file," REMARKS %i patches were applied to the molecule.\n", mol->npatch);
#endif

  }

  ntopo = hasharray_count(mol->defs->topo_hash);
  for ( itopo=0; itopo<ntopo; ++itopo ) {
    topo = &(mol->defs->topo_array[itopo]);
    ntitle_count++;

#if !defined(NOIO)
    fprintf(file," REMARKS topology %s \n", topo->filename);
#endif    

  }

  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    char angles[20], diheds[20];
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    strcpy(angles,"none");
    strcpy(diheds,"");
    if (seg->auto_angles)    strcpy(angles,"angles");
    if (seg->auto_dihedrals) strcpy(diheds,"dihedrals");
    ntitle_count++;
    
#if !defined(NOIO)
    fprintf(file," REMARKS segment %s { first %s; last %s; auto %s %s }\n", seg->segid, seg->pfirst, seg->plast, angles, diheds);
#endif
    
  }

  if (!nopatches) {
  for ( patch = mol->patches; patch; patch = patch->next ) {
    strcpy(defpatch,"");
    if (patch->deflt) strcpy(defpatch,"default");
    npres = patch->npres;
    ipres = 0;
    for ( patchres = patch->patchresids; patchres; patchres = patchres->next ) {
      /* Test the existence of segid:resid for the patch */
      if (!topo_mol_validate_patchres(mol,patch->pname,patchres->segid, patchres->resid)) {
	break;
      }
    }
    if ( patchres ) continue;

    for ( patchres = patch->patchresids; patchres; patchres = patchres->next ) {
      if (ipres==0) {
        ntitle_count++;

#if !defined(NOIO)
        fprintf(file," REMARKS %spatch %s ", defpatch, patch->pname);
#endif      
  
      }
      if (ipres>0 && !ipres%6) {
        ntitle_count++;

#if !defined(NOIO)
        fprintf(file,"\n REMARKS patch ---- ");
#endif

      }

#if !defined(NOIO)
      fprintf(file,"%s:%s  ", patchres->segid, patchres->resid);
      if (ipres==npres-1) fprintf(file,"\n");
#endif

      ipres++;
    }
  }
  }

#if !defined(NOIO)
  fprintf(file,"\n");
  fgetpos(file,&save_pos);
  fsetpos(file,&ntitle_pos);
  fprintf(file,ntitle_fmt,ntitle_count);
  fsetpos(file,&save_pos);
  /* Set 10 Characters for the sections headers for the CHARMM Extended
   * file format. This change is applied for every header
  */
  fprintf(file, charmmext ? "%10d !NATOM\n" : "%8d !NATOM\n",atomid);
#endif

  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      char resid[9];
      res = &(seg->residue_array[ires]);
      strncpy(resid,res->resid,9);
      resid[ charmmext ? 8 : 4 ] = '\0';
      
#if !defined(NEWPSFGEN)  

      if ( charmmfmt ) for ( atom = res->atoms; atom; atom = atom->next ) {
        int idef, typeid;

#else

      if ( charmmfmt ) for ( i = 0; i < res->atomSize; i++ ) {
        int idef, typeid;
        atom = res->atomArray[i];
        atomuint = 0;
#endif

        idef = hasharray_index(mol->defs->type_hash,atom->type);
        if ( idef == HASHARRAY_FAIL ) {
          sprintf(buf,"unknown atom type %s",atom->type);
          print_msg(vdata, v,buf);
          return -3;
        }
        typeid = mol->defs->type_array[idef].id;
        
#if !defined(NOIO)
#if !defined(NEWPSFGEN)
        fprintf(file, ( charmmext ?
                     "%10d %-8s %-8s %-8s %-8s %4d %10.6f    %10.4f  %10d\n" :
                     "%8d %-4s %-4s %-4s %-4s %4d %10.6f    %10.4f  %10d\n" ),
                atom->atomid, seg->segid,resid,res->name,
                atom->name,typeid,atom->charge,atom->mass,0);
                
#else
        if (atom->isdrudlonepair) {
          atomuint = -1;
        } 
          
        if (!numdrude) {
          fprintf(file, ( charmmext ?
                       "%10d %-8s %-8s %-8s %-8s %4d %10.6f    %10.4f  %10d\n" :
                       "%8d %-4s %-4s %-4s %-4s %4d %10.6f    %10.4f  %10d\n" ),
                  atom->atomid, seg->segid,resid,res->name,
                  atom->name,typeid,atom->charge,atom->mass,atomuint);
        } else {
          tmpmass = atom->mass; // the mass of an drude particle host atom was 
                                // already subtracted by 0.4.
          tmpcharge = atom->charge;
          if (!atom->alpha) {
            tmpalpha = 0.0;
            tmpthole = 0.0;
          } else {
            /* Add the drude particles */
            tmpalpha = atom->alpha;
            tmpthole = atom->thole;
          }
          fprintf(file, ( charmmext ?
                       "%10d %-8s %-8s %-8s %-8s %4d %10.6f    %10.4f  %10d %10.6f  %10.6f\n" :
                       "%8d %-4s %-4s %-4s %-4s %4d %10.6f    %10.4f  %10d %10.6f  %10.6f\n" ),
                  atom->atomid, seg->segid, resid, res->name,
                  atom->name, typeid, tmpcharge, tmpmass,atomuint, tmpalpha, tmpthole);
          if (atom->alpha) {
            idef = hasharray_index(mol->defs->type_hash,"DRUD");
            if ( idef == HASHARRAY_FAIL ) {
              sprintf(buf,"unknown atom type %s",atom->type);
              print_msg(vdata, v, buf);
              return -1;
            }
            typeid = mol->defs->type_array[idef].id;
            atomuint = -2;
            fprintf(file, ( charmmext ?
                         "%10d %-8s %-8s %-8s %-8s %4d %10.6f    %10.4f  %10d %10.6f  %10.46f\n" :
                         "%8d %-4s %-4s %-4s %-4s %4d %10.6f    %10.4f  %10d %10.6f  %10.6f\n" ),
                    atom->atomid+1, seg->segid,resid,res->name,
                    atom->dname,typeid,atom->dcharge,0.4,atomuint, 0.0 , 0.0);
          }
        }
#endif
#endif                

                
#if !defined(NEWPSFGEN)  

      } else for ( atom = res->atoms; atom; atom = atom->next ) {

#else

      } else for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        atomuint = 0;
#endif

#if !defined(NOIO)
#if !defined(NEWPSFGEN)
        fprintf(file, ( charmmext ?
                     "%10d %-8s %-8s %-8s %-8s %-6s %10.6f    %10.4f  %10d\n" :
                     "%8d %-4s %-4s %-4s %-4s %-4s %10.6f    %10.4f  %10d\n" ),
                atom->atomid, seg->segid,resid,res->name,
                atom->name,atom->type,atom->charge,atom->mass,0);
                
#else
        /* Add the drude particle */
        if (atom->isdrudlonepair) {
          atomuint = -1;
        } 
        if (!numdrude) {
          fprintf(file, ( charmmext ?
                       "%10d %-8s %-8s %-8s %-8s %-6s %10.6f    %10.4f  %10d\n" :
                       "%8d %-4s %-4s %-4s %-4s %-4s %10.6f    %10.4f  %10d\n" ),
                  atom->atomid, seg->segid,resid,res->name,
                  atom->name,atom->type,atom->charge,atom->mass,atomuint);
        } else {
          tmpmass = atom->mass;
          tmpcharge = atom->charge;
          if (!atom->alpha) {
            tmpalpha = 0.0;
            tmpthole = 0.0;
          } else {
            tmpalpha = atom->alpha;
            tmpthole = atom->thole;
          }
          
          fprintf(file, ( charmmext ?
                        "%10d %-8s %-8s %-8s %-8s %-6s %10.6f    %10.4f  %10d %10.6f  %10.6f\n" :
                        "%8d %-4s %-4s %-4s %-4s %-4s %10.6f    %10.4f  %10d %10.6f  %10.6f\n" ),
                  atom->atomid, seg->segid, resid, res->name,
                  atom->name, atom->type, tmpcharge, tmpmass, atomuint, tmpalpha, tmpthole);
                  
          if (atom->alpha) {
            atomuint = -2;
            fprintf(file, ( charmmext ?
                        "%10d %-8s %-8s %-8s %-8s %-6s %10.6f    %10.4f  %10d %10.6f  %10.6f\n" :
                        "%8d %-4s %-4s %-4s %-4s %-4s %10.6f    %10.4f  %10d %10.6f  %10.6f\n" ),
                    atom->atomid +1, seg->segid,resid,res->name,
                    atom->dname,"DRUD",atom->dcharge,0.4,atomuint, 0.0 , 0.0);
          }
          
        }
#endif

#endif

      }
    }
  }

#if !defined(NOIO)
  fprintf(file,"\n");

  fprintf(file,charmmext ? "%10d !NBOND: bonds\n" :"%8d !NBOND: bonds\n",nbonds);
#endif

  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      
#if !defined(NEWPSFGEN)  

      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( bond = atom->bonds; bond;
                bond = topo_mol_bond_next(bond,atom) ) {
          
#else
      
      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];

        
        
        /* Set the bonds between the lone pairs and drude particles
         * and their hosts if vpbonds = 1
         */
        if (psfcontext->VPBONDS && atom->isdrudlonepair && atom->lonepair->atoms[1]) {
            
#if !defined(NOIO)
          if ( numinline == 4 ) { fprintf(file,"\n");  numinline = 0; }
          fprintf(file, ( charmmext ? " %9d %9d" : " %7d %7d"),
                  atom->atomid,atom->lonepair->atoms[1]->atomid);
#endif
          ++numinline;
          continue;
        } 
        
        for ( bond = atom->bonds; bond;
                bond = topo_mol_bond_next(bond,atom) ) {
          
#endif    
  


#if !defined(NEWPSFGEN)
          if ( bond->atom[0] == atom && ! bond->del ) {
            
            
#if !defined(NOIO)
            if ( numinline == 4 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file, ( charmmext ? " %9d %9d" : " %7d %7d"),
                    atom->atomid,bond->atom[1]->atomid);
#endif
            ++numinline;
          }
#else

          if ( bond->atom[0]->atomid == atom->atomid && ! bond->del ) {
              
            if ( numinline == 4 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file, ( charmmext ? " %9d %9d" : " %7d %7d"),
                    atom->atomid,bond->atom[1]->atomid);
            ++numinline;
            
          }
#endif
        }
        
#if defined(NEWPSFGEN)        
        if (psfcontext->VPBONDS && atom->alpha) {
#if !defined(NOIO)
          if ( numinline == 4 ) { fprintf(file,"\n");  numinline = 0; }
          fprintf(file, ( charmmext ? " %9d %9d" : " %7d %7d"),
                  atom->atomid,atom->atomid +1);
#endif
          ++numinline;
        }
#endif        
        
      }
    }
  }

#if !defined(NOIO)
  fprintf(file,"\n\n");

  fprintf(file,charmmext ? "%10d !NTHETA: angles\n" : "%8d !NTHETA: angles\n",nangls);
#endif

  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      
#if !defined(NEWPSFGEN)

      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( angl = atom->angles; angl;
                angl = topo_mol_angle_next(angl,atom) ) {
          if ( angl->atom[0] == atom && ! angl->del) {

#if !defined(NOIO)
            if ( numinline == 3 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file, ( charmmext ? " %9d %9d %9d" : " %7d %7d %7d"),atom->atomid,
                angl->atom[1]->atomid,angl->atom[2]->atomid);
#endif

            ++numinline;
          }
        }
      }  
#else

      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        for ( angl = atom->angles; angl;
                angl = topo_mol_angle_next(angl,atom) ) {
          if (angl->atom[0]->atomid != atom->atomid || angl->del || 
              angl->atom[1]->isdrudlonepair || 
              angl->atom[2]->isdrudlonepair) continue;
#if !defined(NOIO)
          if ( numinline == 3 ) { fprintf(file,"\n");  numinline = 0; }
          fprintf(file, ( charmmext ? " %9d %9d %9d" : " %7d %7d %7d"),atom->atomid,
              angl->atom[1]->atomid,angl->atom[2]->atomid);
#endif
          ++numinline;
        }
      }
#endif

    }
  }
  
#if !defined(NOIO)
  fprintf(file,"\n\n");

  fprintf(file,charmmext ? "%10d !NPHI: dihedrals\n" : "%8d !NPHI: dihedrals\n",ndihes);
#endif

  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      
#if !defined(NEWPSFGEN)

      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( dihe = atom->dihedrals; dihe;
                dihe = topo_mol_dihedral_next(dihe,atom) ) {
          if ( dihe->atom[0] == atom && ! dihe->del ) {

#if !defined(NOIO)
            if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file, ( charmmext ? " %9d %9d %9d %9d" : " %7d %7d %7d %7d"),atom->atomid,
                dihe->atom[1]->atomid,dihe->atom[2]->atomid,
                dihe->atom[3]->atomid);
#endif
            ++numinline;
          }
        }
#else

      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        if (atom->isdrudlonepair) continue;
        for ( dihe = atom->dihedrals; dihe;
                dihe = dihe->next ) {
          if (dihe->del || dihe->atom[1]->isdrudlonepair || 
              dihe->atom[2]->isdrudlonepair || dihe->atom[0]->del ||
              dihe->atom[1]->del || dihe->atom[2]->del) 

            continue;
#if !defined(NOIO)
          if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
          fprintf(file, ( charmmext ? " %9d %9d %9d %9d" : " %7d %7d %7d %7d"),atom->atomid,
              dihe->atom[0]->atomid,dihe->atom[1]->atomid,
              dihe->atom[2]->atomid);
#endif
          ++numinline;
        }
#endif

      }
    }
  }

#if !defined(NOIO)
  fprintf(file,"\n\n");

  fprintf(file,charmmext ? "%10d !NIMPHI: impropers\n" : "%8d !NIMPHI: impropers\n",nimprs);
#endif

  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)

      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( impr = atom->impropers; impr;
                impr = topo_mol_improper_next(impr,atom) ) {
          if ( impr->atom[0] == atom && ! impr->del ) {
#if !defined(NOIO)
            if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file, ( charmmext ? " %9d %9d %9d %9d" : " %7d %7d %7d %7d"),atom->atomid,
                impr->atom[1]->atomid,impr->atom[2]->atomid,
                impr->atom[3]->atomid);
#endif
          ++numinline;
          }
#else

      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        if (atom->isdrudlonepair) continue;
        for ( impr = atom->impropers; impr;
                impr = impr->next) {
          if (impr->del || impr->atom[1]->isdrudlonepair || 
              impr->atom[2]->isdrudlonepair || impr->atom[0]->del || 
              impr->atom[1]->del || impr->atom[2]->del) 

            continue;
#if !defined(NOIO)
            if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file, ( charmmext ? " %9d %9d %9d %9d" : " %7d %7d %7d %7d"),atom->atomid,
                impr->atom[0]->atomid,impr->atom[1]->atomid,
                impr->atom[2]->atomid);
#endif
          ++numinline;
#endif
        }
      }
    }
  }

#if !defined(NOIO)
  fprintf(file,"\n\n");

  fprintf(file,charmmext ? "%10d !NDON: donors\n\n\n" : "%8d !NDON: donors\n\n\n",0);
  fprintf(file,charmmext ? "%10d !NACC: acceptors\n\n\n" : "%8d !NACC: acceptors\n\n\n",0);
  fprintf(file,charmmext ? "%10d !NNB\n" : "%8d !NNB\n",nexcls);
#endif

  /* Print atom numbers for exclusions */
  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)

      for ( atom = res->atoms; atom; atom = atom->next ) {
        
#else

      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        //repeat exclusions for drude particles to match parent
        for ( j = 0; j < (atom->alpha ? 2 : 1); j++ ) {
#endif
        
          for ( excl = atom->exclusions; excl;
              excl = topo_mol_exclusion_next(excl,atom) ) {
            if ( excl->atom[0]->atomid == atom->atomid && ! excl->del) {

#if !defined(NOIO)
              if ( numinline == 8 ) { fprintf(file,"\n");  numinline = 0; }
              fprintf(file,(charmmext?" %9d":" %7d"),excl->atom[1]->atomid);
#endif

              ++numinline;
            
#if defined(NEWPSFGEN)
            }
#endif
          }
        }
      }
    }
  }

#if !defined(NOIO)
  fprintf(file,"\n");
#endif

  /* Print exclusion indices for every atom */
  nexcls = 0;
  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)

      for ( atom = res->atoms; atom; atom = atom->next ) {
        
#else

      for ( i = 0; i < res->atomSize; i++ ) {
        atom = res->atomArray[i];
        //repeat exclusions for drude particles to match parent
        for ( j = 0; j < (atom->alpha ? 2 : 1); j++ ) {
#endif

          for ( excl = atom->exclusions; excl;
                  excl = topo_mol_exclusion_next(excl,atom) ) {
            if ( excl->atom[0] == atom && ! excl->del ) {
                ++nexcls;
            }
          }
        
#if !defined(NOIO)
          if ( numinline == 8 ) { fprintf(file,"\n");  numinline = 0; }
          fprintf(file,(charmmext?" %9d":" %7d"),nexcls);
#endif


        ++numinline;
#if defined(NEWPSFGEN)
        }
#endif
      }
    }
  }

#if !defined(NOIO)
  fprintf(file,"\n\n");

  fprintf(file,(charmmext?"%10d %9d !NGRP\n%10d%10d%10d\n\n":"%8d %7d !NGRP\n%8d%8d%8d\n\n"),1,0,0,0,0);
#endif

#if defined(NEWPSFGEN)
  /*By Brian Radak
  * Requirements for NAMD as I understand them:
  * 1) LP has 0 mass and follows parent atoms (the "LP hosts", in this case CL and C6)
  *
  * 2) bond between LP and parent should be optional in the bonds section, but 
  *    this might be necessary * for building migration groups properly?
  *
  * 3) no angles, dihedrals, etc. should contain LP, this might just be because those 
  *    types are not * defined

  * 4) the new "!NUMLP NUMLPH" section in the PSF should exist. The presence of 
  *    this section * automatically toggles the now deprecated "lonepairs on" 
  *    keyword and is also required for "drude * on". Format is as follows:
  * =========
  * <# of lonepairs> <# of lonepairs + # of lphosts> !NUMLP NUMLPH
  * .
  * :
  * <# of lphosts for this lonepair> <the pointer for this lonepair> F <distance> <angle> <dihedral>
  * .
  * :
  * <lonepair indices>
  * <lonepair index> <lphost1> <lphost2> .... <lonepair index> ...
  * ==========

  * Note the mislabeling of NUMLPH in the comment, the actual number of hosts is 
  * the second number * minus the first.

  * For each entry, a collinear lonepair has 2 hosts and a bisectory has 3. I 
  * forget what F mean * ("fixed"?) this is the only option I've seen and the 
  * only one NAMD accepts. All entries require a * <distance>, <angle>, <dihedral>  
  * specification. For collinear lonepairs dihedral is read, but * ignored. <angle> 
  * is interpreted as a scale parameter that shifts the origin for the colinear x distance 
  * back along the bond between the parent atom and the other host. A value of 0.0 means 
  * that * the origin of the distance is the parent atom.
  *
  * The lonepairs are essentially double indexed (starting at one in Fortran style), 
  * once as an atom * and again within the list of lonepairs and hosts. So for the example, 
  * the first lonepair has index * 1, which refers to the 13, which is the lonepair 
  * atom index. The 2 for that entry indicates two * extra host entries (12 and 11, 
  * lonepair indices 2 and 3) and because it is collinear, the lonepair * is bound to atom 12.

  * An additional lonepair would have index 4, add either 2 or 3 more hosts to the 
  * index list and thus * add 3 or 4 more entries. The index list should wrap to the 
  * next line after every 8 entries.
  * Example, lets add two 5 point SWM4 waters with bisector lonepairs after the chorobenzene. 
  * The new section would be:

  * 1         3 !NUMLP NUMLPH
  * 2         1   F   1.64000       0.00000       0.00000    
  * 3         4   F  -0.240345     0.00000       0.00000
  * 3         8   F  -0.240345     0.00000       0.00000
  * 13        12        11        16       14      17        18        21
  * 19        22        23
  */
  if (numlonepairs) {
    
#if !defined(NOIO)
  fprintf(file,charmmext ? "%10d %9d !NUMLP NUMLPH\n" : "%8d %7d !NUMLP NUMLPH\n", numlonepairs,numlphosts + numlonepairs);
#endif

    numlonepairs = 1;
    for ( iseg=0; iseg<nseg; ++iseg ) {
      seg = mol->segment_array[iseg];
      if (! seg) continue;
      nres = hasharray_count(seg->residue_hash);
      for ( ires=0; ires<nres; ++ires ) {
        res = &(seg->residue_array[ires]);
        for ( i = 0; i < res->atomSize; i++ ) {
          atom = res->atomArray[i];
          if (atom->isdrudlonepair) {
            switch (atom->lonepair->lptype) {
              case COLINEARLP: 
                numlphosts = 2;
                lpunit = 1.0;
                break;
              case RELATIVELP:
              case BISECTORLP:
                numlphosts = 3;
                lpunit = 180.0/M_PI;
                break;
              default:
                numlphosts = 3;
                lpunit = 180.0/M_PI;
                break;
            }
#if !defined(NOIO)
            fprintf(file, ("%8d %8d %-6s %10.4f %10.4f %10.4f\n"),
                          numlphosts,  numlonepairs,"F", atom->lonepair->distance,
                        atom->lonepair->angle * lpunit, atom->lonepair->dihedral * lpunit);
#endif
            numlonepairs += numlphosts + 1;
          }
        }
      }
    }
    numinline = 0;
    numlphosts = 0;
    for ( iseg=0; iseg<nseg; ++iseg ) {
      seg = mol->segment_array[iseg];
      if (! seg) continue;
      nres = hasharray_count(seg->residue_hash);
      for ( ires=0; ires<nres; ++ires ) {
        res = &(seg->residue_array[ires]);
        for ( i = 0; i < res->atomSize; i++ ) {
          atom = res->atomArray[i];
          if (atom->isdrudlonepair) {
            switch (atom->lonepair->lptype) {
              case COLINEARLP: 
                numlphosts = 3;
                break;
              case RELATIVELP:
              case BISECTORLP:
                numlphosts = 4;
                break;
              default:
                numlphosts = 4;
                break;
            }
            for (j = 0 ; j < numlphosts; j++) {
              
#if !defined(NOIO)
              fprintf(file, ( charmmext ? " %9d" : " %7d"),
              atom->lonepair->atoms[j]->atomid);
#endif

              ++numinline;
              if ( numinline == 8 ) { fprintf(file,"\n");  numinline = 0; }
            }
          }
        }
      }
    }
    fprintf(file,"\n\n");
  }
  
  
  if (numaniso || numdrude) {
#if !defined(NOIO)
    fprintf(file,charmmext ? "%10d !NUMANISO\n" :"%8d !NUMANISO\n",numaniso);
#endif
    if (numaniso) {
      for ( iseg=0; iseg<nseg; ++iseg ) {
        seg = mol->segment_array[iseg];
        if (! seg) continue;
        nres = hasharray_count(seg->residue_hash);
        for ( ires=0; ires<nres; ++ires ) {
          res = &(seg->residue_array[ires]);
          if (res->numaniso ) {
            for ( aniso = res->aniso; aniso; aniso = aniso->next ) { 
              if (aniso->del) continue;
            
#if !defined(NOIO)
              fprintf(file, "%20.4f %14.4f %14.4f\n",
                            aniso->k11,  aniso->k22, aniso->k33);
#endif

            }
          }
        }
      }
      numinline = 0;
      for ( iseg=0; iseg<nseg; ++iseg ) {
        seg = mol->segment_array[iseg];
        if (! seg) continue;
        nres = hasharray_count(seg->residue_hash);
        for ( ires=0; ires<nres; ++ires ) {
          res = &(seg->residue_array[ires]);
          for ( aniso = res->aniso; aniso; aniso = aniso->next ) { 
            if (aniso->del) continue;
#if !defined(NOIO)
            fprintf(file, ( charmmext ? " %9d %9d %9d %9d" : " %7d %7d %7d %7d"),
            aniso->atoms[0]->atomid,  aniso->atoms[1]->atomid, 
            aniso->atoms[2]->atomid, aniso->atoms[3]->atomid);
#endif
            ++numinline;
            if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
          }
        }
      }
    }
#if !defined(NOIO)
    fprintf(file,"\n\n");
#endif

  }
  
  
#endif
  if ( ! nocmap ) {
    
#if !defined(NOIO)
    fprintf(file,charmmext ? "%10d !NCRTERM: cross-terms\n":"%8d !NCRTERM: cross-terms\n",ncmaps);
#endif
    for ( iseg=0; iseg<nseg; ++iseg ) {
      seg = mol->segment_array[iseg];
      if (! seg) continue;
      nres = hasharray_count(seg->residue_hash);
      for ( ires=0; ires<nres; ++ires ) {
        res = &(seg->residue_array[ires]);

#if !defined(NEWPSFGEN)

        for ( atom = res->atoms; atom; atom = atom->next ) {
        
#else

        for ( i = 0; i < res->atomSize; i++ ) {
          atom = res->atomArray[i];
        
#endif

          for ( cmap = atom->cmaps; cmap;
                  cmap = topo_mol_cmap_next(cmap,atom) ) {
            if ( cmap->atom[0] == atom && ! cmap->del ) {

#if !defined(NOIO)
              fprintf(file,( charmmext ? " %9d %9d %9d %9d %9d %9d %9d %9d\n"
                         : " %7d %7d %7d %7d %7d %7d %7d %7d\n"),atom->atomid,
                  cmap->atom[1]->atomid,cmap->atom[2]->atomid,
                  cmap->atom[3]->atomid,cmap->atom[4]->atomid,
                  cmap->atom[5]->atomid,cmap->atom[6]->atomid,
                  cmap->atom[7]->atomid);
#endif                  

            }
          }
        }
      }
    }

#if !defined(NOIO)
    fprintf(file,"\n");
#endif

  }

  return 0;
}




