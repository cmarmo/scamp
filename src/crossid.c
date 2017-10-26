/*
 *    crossid.c
 *
 * Manage source cross-identifications.
 *
 * This file part of: SCAMP
 *
 * Copyright:  (C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
 *
 * License:  GNU General Public License
 *
 * SCAMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * SCAMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 * Last modified:  13/09/2012
 *
 */
#include "define.h"
#include "globals.h"
#include "crossid.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "match.h"
#include "misc.h"
#include "prefs.h"
#include "samples.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/*
 * Store crossid calculus states
 */
struct cross_state {
    double lngmin;
    double lngmax;
    double latmin;
    double latmax;
    double projmin[NAXIS];
    double projmax[NAXIS];
};

/* 
 * Compute the largest possible error in pixels allowed in previous matching 
 */
double 
__get_limit(fgroupstruct *fgroup, 
            double        tolerance) 
{
    double lim_tmp;
    double lim_result = 0.0;
    int i;

    for (i=0; i < fgroup->naxis; i++) {
        lim_tmp = tolerance / fgroup->meanwcsscale[i];
        if (lim_tmp > lim_result)
            lim_result = lim_tmp;
    }

    return lim_result;
}


/* 
 * Sort samples to accelerate further processing and reset pointers 
 */
void 
__sort_samples(fgroupstruct *fgroup) 
{
    fieldstruct **field = fgroup->field;
    setstruct    *set;
    int i, j;

    for (i=0; i < fgroup->nfield; i++) {

        field[i]->prevfield = NULL;
        field[i]->nextfield = NULL;

        for (j=0; j < field[i]->nset; j++) {

            set = field[i]->set[j];

            sort_samples(set);
            unlink_samples(set);

        }
    }
}

/*
 * Exclude non overlapping frames.
 * Returns 1 if it overlaps, -1 if not.
 */
int
__cross_overlaps(setstruct *set_a, 
                setstruct *set_b,
                int        lng, 
                int        lat, 
                int        naxis,
                double     limit,
                struct cross_state *state)
{
    int i;

    if (lng != lat) {

        if ( set_a->projposmin[lng] > (state->lngmax = set_b->projposmax[lng] + limit) ||
            (state->lngmin = set_b->projposmin[lng] - limit) > set_a->projposmax[lng]  ||
             set_a->projposmin[lat] > (state->latmax = set_b->projposmax[lat] + limit) ||
            (state->latmin = set_b->projposmin[lat] - limit)> set_a->projposmax[lat]) {
            return -1;
        }

    } else {

        for (i=0; i<naxis; i++) {
            if ( set_a->projposmin[i] > (state->projmax[i]=set_b->projposmax[i]+limit) ||
                (state->projmin[i] = set_b->projposmin[i]-limit)> set_a->projposmax[i]) 
            {   
                return -1;
            }   
        }   
    }

    return 1;
}

/*
 * Crossid algorythm.
 */
void
__cross_set(setstruct *set_a, 
            setstruct *set_b,
            wcsstruct *wcs,
            double     limit,
            int        naxis,
            int        lng,
            int        lat)
{
    int status;
    struct cross_state state;
    state.lngmin = 0.0;
    state.lngmax = 0.0;
    state.latmin = 0.0;
    state.latmax = 0.0;

    status = __cross_overlaps(set_a, set_b, lng, lat, naxis, limit, &state);
    if (status < 0)
        return;
    
    // STOP HERE, restart at "for(nsamp-set1->nsample;nsam--;samp1++..."
}

/*
 * Loop over fields and apply crossid algorythm
 */
void
__loop_crossid(fgroupstruct *fgroup, double limit)
{
    int i, j, k, l;

    struct set *set_a;
    struct set *set_b;

    fieldstruct **fields = fgroup->field;
    wcsstruct    *wcs    = fgroup->wcs;
    int naxis            = fgroup->naxis;
    int lng              = fgroup->lng;
    int lat              = fgroup->lat;

    /*
     * With this, we are crossing every sets together once.
     * Can this be parallelized?
     */
    for (i=0; i < fgroup->nfield; i++) {
        for (j=0; j < fields[i]->nset; j++) {
            for (k=0; k < i; k++) {
                for (l=0; l < fields[k]->nset; l++) {

                    set_a = fields[i]->set[j];
                    set_b = fields[k]->set[l];

                    __cross_set(set_a, set_b, wcs, limit, naxis, lng, lat);

                }
            }
        }
    }
}


/**
 * @fn void crossid_fgroup(fgroupstruct *fgroup, fieldstruct *reffield, double tolerance)
 * @brief Perform source cross-identifications in a group of fields.
 * @param fgroup ptr to the group of fields,
 * @param reffield ptr to the reference field,
 * @param tolerance Tolerance (in deg. if angular coordinates).
 * @remarks Uses the global preferences.
 * @author E. Bertin (IAP)
 * @version 13/09/2012
 */
void crossid_fgroup(fgroupstruct *fgroup, 
                    fieldstruct  *reffield,
                    double        tolerance)
{
    fieldstruct **field, *field1, *field2;
    wcsstruct *wcs;
    setstruct **pset1,**pset2, **pset,
              *set1, *set2, *set;
    samplestruct *samp1, *nexsamp1, *samp2, *samp2b, *samp2min, 
                 *prevsamp2, *nextsamp1, *nextsamp2;
    double projmin2[NAXIS], projmax2[NAXIS], *proj1,
           lng1,lat1, latmin1, latmax1, lngmin2, lngmax2, latmin2, latmax2,
           dlng, dlat, dx, rlimit, r2, r2n, r2p, r2min;
    float fmax;
    int i, f, f1, f2, s1, s2, nset1, nset2, nsamp, nsamp2, nsamp2b,
        s, nfield, naxis, lng,lat, yaxis;

    lng1 = lngmin2 = lngmax2 = latmin2 = latmax2 = 0.0;
    field1 = NULL; /* to avoid gcc -Wall warnings */
    proj1  = NULL; /* to avoid gcc -Wall warnings */
    field   = fgroup->field;
    nfield  = fgroup->nfield;
    naxis   = fgroup->naxis;
    lng     = fgroup->lng;
    lat     = fgroup->lat;
    wcs     = fgroup->wcs;


    rlimit = __get_limit(fgroup, tolerance);
    __sort_samples(fgroup);

    /*
     * Now start the real cross-id loop 
     */

    __loop_crossid(fgroup, rlimit);

    /* Foreach fealds */
    for (f1=1; f1<nfield; f1++)
    {
        field1 = field[f1];
        pset1  = field1->set;
        set1   = *(pset1++);
        nset1  = field1->nset;

        /* Foreach sets */
        for (s1=nset1; s1--; set1=*(pset1++))
        {

            /* Foreach precedent fields */
            for (f2=0; f2<f1; f2++)
            {
                field2  = field[f2];
                pset2   = field2->set;
                set2    = *(pset2++);
                nset2   = field2->nset;

                /* Foreach sets of the precedent field */
                for (s2=nset2; s2--; set2=*(pset2++))
                {
                    /*-------- Exclude non-overlapping frames */
                    if (lng != lat)
                    {
                        if (
                                set1->projposmin[lng] > (lngmax2=set2->projposmax[lng] + rlimit) ||
                                (lngmin2=set2->projposmin[lng] - rlimit) > set1->projposmax[lng] ||
                                set1->projposmin[lat] > (latmax2=set2->projposmax[lat] + rlimit) ||
                                (latmin2=set2->projposmin[lat] - rlimit)> set1->projposmax[lat]) 
                        {
                            continue;
                        }
                    } 
                    else 
                    {
                        for (i=0; i<naxis; i++) {
                            if (
                                    set1->projposmin[i] > (projmax2[i]=set2->projposmax[i]+rlimit) ||
                                    (projmin2[i]=set2->projposmin[i]-rlimit)> set1->projposmax[i]) 
                            {
                                continue;
                            }
                        }
                    }
                    samp1   = set1->sample;
                    samp2b  = set2->sample;
                    nsamp2b = set2->nsample;

                    /* Foreach ??? */
                    for (nsamp=set1->nsample; nsamp--; samp1++)
                    {
                        if (lat!=lng)
                        {
                            lng1  = samp1->projpos[lng];
                            lat1  = samp1->projpos[lat];
                            yaxis = lat;
                            /*------------ Jump over sources in the non-overlapping region */
                            if (lat1<latmin2 || lat1>latmax2 || lng1<lngmin2 || lng1>lngmax2)
                                continue;
                        }
                        else
                        {
                            proj1 = samp1->projpos;
                            for (i=0; i<naxis; i++) {
                                if (proj1[i] < projmin2[i] || proj1[i]>projmax2[i])
                                    continue;
                            }
                            lat1 = (naxis<2) ? proj1[yaxis=0] : proj1[yaxis=1];
                        }

                        latmin1  = lat1-rlimit;
                        latmax1  = lat1+rlimit;
                        r2min    = rlimit*rlimit;
                        samp2min = NULL;
                        samp2    = samp2b;
                        /*---------- Jump over sources that can't match in y */

                        /* Foreach What??? */
                        for (nsamp2=nsamp2b; nsamp2-- && samp2->projpos[yaxis] < latmin1; samp2++);

                        samp2b  = samp2;
                        nsamp2b = ++nsamp2;
                        for (; nsamp2-- && samp2->projpos[yaxis]<latmax1; samp2++)
                        {
                            if (samp2->nextsamp && samp2->nextsamp->set->field!=field1)
                                continue;
                            if (lat!=lng)
                            {
                                dlng = lng1 - samp2->projpos[lng];
                                dlat = lat1 - samp2->projpos[lat];
                                r2   = dlng*dlng + dlat*dlat;
                            }
                            else
                            {
                                r2 = 0.0;
                                for (i=0; i<naxis; i++)
                                {
                                    dx  = proj1[i] - samp2->projpos[i];
                                    r2 += dx*dx;
                                }
                            }
                            /*------------ Finally select the closest source within the search disk */
                            if (r2<r2min)
                            {
                                r2min    = r2;
                                samp2min = samp2;
                            }
                        }
                        if (samp2min)
                        {
                            r2p = r2n = BIG;
                            if ((prevsamp2=samp1->prevsamp))
                            {
                                /*-------------- Check if it is a better match than the previous one */
                                if (lat!=lng)
                                {
                                    dlng = prevsamp2->projpos[lng] - samp1->projpos[lng];
                                    dlat = prevsamp2->projpos[lat] - samp1->projpos[lat];
                                    r2p  = dlng*dlng + dlat*dlat;
                                }
                                else
                                {
                                    r2p = 0.0;
                                    for (i=0; i<naxis; i++)
                                    {
                                        dx   = prevsamp2->projpos[i] - samp1->projpos[i];
                                        r2p += dx*dx;
                                    }
                                }
                            }
                            if ((nextsamp1=samp2min->nextsamp))
                            {
                                /*-------------- Check if it is a better match than the previous one */
                                if (lat!=lng)
                                {
                                    dlng = samp2min->projpos[lng] - nextsamp1->projpos[lng];
                                    dlat = samp2min->projpos[lat] - nextsamp1->projpos[lat];
                                    r2n  = dlng*dlng + dlat*dlat;
                                }
                                else
                                {
                                    r2n = 0.0;
                                    for (i=0; i<naxis; i++)
                                    {
                                        dx   = samp2min->projpos[i] - nextsamp1->projpos[i];
                                        r2n += dx*dx;
                                    }
                                }
                            }
                            if (r2min<r2p && r2min<r2n)
                                /*------------ unlink from previous match if this is a better match */
                            {
                                if (prevsamp2)
                                    prevsamp2->nextsamp = NULL;

                                if (nextsamp1)
                                    nextsamp1->prevsamp = NULL;

                                samp1->prevsamp    = samp2min;
                                samp2min->nextsamp = samp1;

                            }
                        }
                    }
                }
            }
        }
    }

    /* Now bring also the reference field samples to the common projection */
    /* Sort samples to accelerate further processing and reset pointers */
    if (reffield)
    {
        set1 = reffield->set[0];
        sort_samples(set1);
        unlink_samples(set1);

        for (f2=0; f2<nfield; f2++)
        {
            field2 = field[f2];
            pset2  = field2->set;
            set2   = *(pset2++);
            nset2  = field2->nset;

            for (s2=nset2; s2--; set2=*(pset2++))
            {
                /*---------- Exclude non-overlapping frames */
                if (lng != lat)
                {
                    if (set1->projposmin[lng] > (lngmax2=set2->projposmax[lng]+rlimit)
                            || (lngmin2=set2->projposmin[lng]-rlimit) > set1->projposmax[lng]
                            || set1->projposmin[lat] > (latmax2=set2->projposmax[lat]+rlimit)
                            || (latmin2=set2->projposmin[lat]-rlimit)> set1->projposmax[lat])
                        continue;
                }
                else
                {
                    for (i=0; i<naxis; i++)
                    {
                        if (set1->projposmin[i] > (projmax2[i]=set2->projposmax[i]+rlimit)
                                || (projmin2[i]=set2->projposmin[i]-rlimit) >set1->projposmax[i])
                            continue;
                    }
                }
                samp1   = set1->sample;
                samp2b  = set2->sample;
                nsamp2b = set2->nsample;
                for (nsamp=set1->nsample; nsamp--; samp1++)
                {
                    if (lat!=lng)
                    {
                        lng1  = samp1->projpos[lng];
                        lat1  = samp1->projpos[lat];
                        yaxis = lat;
                        /*---------- Jump over sources in the non-overlapping region */
                        if (lat1<latmin2 || lat1>latmax2 || lng1<lngmin2 || lng1>lngmax2)
                            continue;
                    }
                    else
                    {
                        proj1 = samp1->projpos;
                        for (i=0; i<naxis; i++)
                            if (proj1[i] < projmin2[i] || proj1[i]>projmax2[i])
                                continue;
                        lat1 = (naxis<2) ? proj1[yaxis=0] : proj1[yaxis=1];
                    }
                    latmin1  = lat1-rlimit;
                    latmax1  = lat1+rlimit;
                    r2min    = rlimit*rlimit;
                    fmax     = 0.0;
                    samp2min = NULL;
                    samp2    = samp2b;
                    /*-------- Jump over sources that can't match in y */
                    for (nsamp2=nsamp2b; nsamp2-- && samp2->projpos[yaxis]<latmin1; samp2++);

                    samp2b  = samp2;
                    nsamp2b = ++nsamp2;
                    for (; nsamp2-- && samp2->projpos[yaxis]<latmax1; samp2++)
                    {
                        if (samp2->prevsamp)
                            continue;
                        if (lat!=lng)
                        {
                            dlng = lng1 - samp2->projpos[lng];
                            dlat = lat1 - samp2->projpos[lat];
                            r2   = dlng*dlng + dlat*dlat;
                        }
                        else
                        {
                            r2 = 0.0;
                            for (i=0; i<naxis; i++)
                            {
                                dx  = proj1[i] - samp2->projpos[i];
                                r2 += dx*dx;
                            }
                        }
                        /*---------- Finally select the closest source within the search disk */
                        if (r2<r2min)
                        {
                            r2min    = r2;
                            samp2min = samp2;
                        }
                    }
                    if (samp2min)
                    {
                        r2n = BIG;
                        if ((nextsamp2=samp1->nextsamp))
                        {
                            /*------------ Check if it is a better match than the previous one */
                            if (lat!=lng)
                            {
                                dlng = nextsamp2->projpos[lng] - samp1->projpos[lng];
                                dlat = nextsamp2->projpos[lat] - samp1->projpos[lat];
                                r2n  = dlng*dlng + dlat*dlat;
                            }
                            else
                            {
                                r2n = 0.0;
                                for (i=0; i<naxis; i++)
                                {
                                    dx   = nextsamp2->projpos[i] - samp1->projpos[i];
                                    r2n += dx*dx;
                                }
                            }
                        }
                        if (r2min<r2n)
                            /*---------- unlink from previous match if this is a better match */
                        {
                            if (nextsamp2)
                                nextsamp2->prevsamp = NULL;
                            samp1->nextsamp    = samp2min;
                            samp2min->prevsamp = samp1;
                        }
                    }
                }
            }
        }
    }

    return;
}


/**
 * @fn void recenter_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
 * @brief Perform field recentering with respect to a reference catalog in a
 *  group of fields.
 * @param fgroup ptr to the group of fields,
 * @param reffield ptr to the reference field.
 * @returns void
 * @remarks Uses the global preferences.
 * @author E. Bertin (IAP)
 * @version 09/06/2011
 */
void recenter_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
{
    fieldstruct *field;
    setstruct **sets, *set;
    samplestruct *samp, *samp2;
    double *offsetbuf[NAXIS],
           offset[NAXIS], rawpos[NAXIS], wcspos[NAXIS], dwcspos[NAXIS];  
    int  d,f,s, naxis, nsamp, o,omax;

    NFPRINTF(OUTPUT, "Re-centering fields...");

    set   = reffield->set[0];
    naxis = fgroup->naxis;
    omax  = 0;
    for (f=0; f<fgroup->nfield; f++)
    {
        o = 0;
        field = fgroup->field[f];
        samp = set->sample;
        for (nsamp=set->nsample; nsamp--; samp++)
        {
            for (samp2=samp; (samp2=samp2->nextsamp);)
            {
                if (samp2->set && samp2->set->field == field)
                {
                    if (o>=omax)
                    {
                        omax += 8192;
                        if (o)
                            for (d=0; d<naxis; d++)
                            {
                                QREALLOC(offsetbuf[d], double, omax);
                            }
                        else
                            for (d=0; d<naxis; d++)
                            {
                                QMALLOC(offsetbuf[d], double, omax);
                            }
                    }
                    for (d=0; d<naxis; d++)
                        offsetbuf[d][o] = samp2->projpos[d] - samp->projpos[d];
                    o++;
                }
            }
        }
        /*-- Compute the median reprojected shift in each dimension */
        for (d=0; d<naxis; d++)
            offset[d] = fast_median(offsetbuf[d], o);
        /*-- Convert it to a shift in world coordinates */
        for (d=0; d<naxis; d++)
            rawpos[d] = field->set[0]->wcs->crpix[d] - offset[d];
        raw_to_wcs(field->set[0]->wcs, rawpos, wcspos);
        for (d=0; d<naxis; d++)
            dwcspos[d] = wcspos[d] - field->set[0]->wcs->crval[d];
        sets = field->set;
        for (s=0; s<field->nset; s++)
            update_wcsll(sets[s]->wcs, dwcspos[set->lng], dwcspos[set->lat]);
    }

    for (d=0; d<naxis; d++)
        free(offsetbuf[d]);

    return;
}


/**
 * @fn int check_fieldoverlap(fieldstruct *field1, fieldstruct *field2)
 * @brief Check if two fields overlap or not.
 * @param field1 ptr to the first field,
 * @param field2 ptr to the second field.
 * @return integer 1 if they overlap, 0 otherwise.
 * @remarks
 * @author E. Bertin (IAP)
 * @version 07/02/2005
 */
int check_fieldoverlap(fieldstruct *field1, fieldstruct *field2)

{
    setstruct **pset,
              *set;
    samplestruct *samp,*samp2;
    int  n,s;

    pset = field1->set;
    set  = *(pset++);
    for (s=field1->nset; s--; set=*(pset++))
    {
        samp = set->sample;
        for (n=set->nsample; n--; samp++)
        {
            if (samp->nextsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->nextsamp)) {
                    if (samp2->set->field == field2)
                        return 1;
                }
            }
            if (samp->prevsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->prevsamp)) {
                    if (samp2->set->field == field2)
                        return 1;
                }
            }
        }
    }

    /* No link found between both fields */
    return 0;
}


/**
 * @fn int check_fieldphotomoverlap(fieldstruct *field, int instru)
 * @brief Check if a field overlaps a photometric field or not.
 * @param field ptr to the field to check,
 * @param instru photometric instrument index.
 * @returns Photometric code (1 for genuine, 2 for dummy) if it overlaps, 0
 * otherwise.
 * @remarks
 * @author E. Bertin (IAP)
 * @version 25/02/2005
 */
int check_fieldphotomoverlap(fieldstruct *field, int instru)
{
    setstruct **pset,
              *set;
    samplestruct *samp,*samp2;
    int  n,s;

    pset = field->set;
    set  = *(pset++);
    for (s=field->nset; s--; set=*(pset++))
    {
        samp = set->sample;
        for (n=set->nsample; n--; samp++)
        {
            if (samp->nextsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->nextsamp)) {
                    if (samp2->set->field->photomflag
                            && samp2->set->field->photomlabel==instru)
                        return samp2->set->field->photomflag;
                }
            }
            if (samp->prevsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->prevsamp)) {
                    if (samp2->set->field->photomflag
                            && samp2->set->field->photomlabel==instru)
                        return samp2->set->field->photomflag;
                }
            }
        }
    }

    /* No photometric field found */
    return 0;
}


