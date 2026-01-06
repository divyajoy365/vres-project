/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/


#include "./custom.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


void contact_function( Cell* pMe, Phenotype& phenoMe, Cell* pOther, Phenotype& phenoOther, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );
void tcell_update_velocity( Cell* pCell, Phenotype& phenotype, double dt );
void tcell_division_function( Cell* pParent, Cell* pDaughter );

std::vector<std::string> coloring_function( Cell* pCell );

void create_cell_types( void )
{
    if( parameters.ints.find_index("random_seed") != -1 )
    {
        SeedRandom( parameters.ints("random_seed") );
    }

    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity       = standard_update_cell_velocity;

    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype      = NULL;
    cell_defaults.functions.custom_cell_rule      = custom_function;

    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane          = NULL;

    cell_defaults.functions.contact_function = NULL;

    if( parameters.doubles.find_index("attachment_elastic_constant") != -1 )
    {
        cell_defaults.phenotype.mechanics.attachment_elastic_constant =
            parameters.doubles("attachment_elastic_constant");
    }

 
    initialize_cell_definitions_from_pugixml();

    // T cell
    Cell_Definition* pT = find_cell_definition( "T cell" );
    if( pT == nullptr )
    {
        std::cout << "Warning: no cell definition named 'T cell' was found.\n";
    }
    else
    {
        pT->functions.update_velocity = tcell_update_velocity;
        pT->functions.contact_function = NULL;
        pT->functions.cell_division_function  = tcell_division_function;
    }

    // rod cell
    Cell_Definition* pRod = find_cell_definition( "rod" );
    if( pRod != nullptr )
    {

        pRod->phenotype.motility.is_motile = false;
        int idx = pRod->custom_data.find_variable_index( "max_attachments" );
        if( idx < 0 )
        {
            pRod->custom_data.add_variable( "max_attachments", "dimensionless", 2.0 );
        }
        else
        {
            pRod->custom_data[ "max_attachments" ] = 2.0;
        }

        pRod->functions.contact_function = contact_function;
    }
    else
    {
        std::cout << "no rod definition";
    }

    build_cell_definitions_maps();
    setup_signal_behavior_dictionaries();
    setup_cell_rules();
    display_cell_definitions( std::cout );
    return;
}

void setup_microenvironment( void )
{
    initialize_microenvironment();
    return;
}

void setup_tissue( void )
{

    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];


    Zmin = 0.0;
    Zmax = 0.0;

    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;


    load_cells_from_pugixml();

    Cell_Definition* pT = find_cell_definition( "T cell" );

    if( pT != nullptr )
    {
        int n_T = 0;
        if( parameters.ints.find_index("number_of_cells") != -1 )
        {
            n_T = parameters.ints("number_of_cells");
        }

        for( int n = 0; n < n_T; n++ )
        {
            // asssign random positions 
            Cell* pC = create_cell( *pT );
            double x = Xmin + UniformRandom() * Xrange;
            
            double y = Ymin + UniformRandom() * Yrange;
            double z = 0.0;

            pC->assign_position( x, y, z );
        }
    }
    else
    {
        std::cout << "setup_tissue: no t cell.\n";
    }

    // create rod cells
    Cell_Definition* pRodDef = find_cell_definition( "rod" );
    if( pRodDef == nullptr )
    {
        std::cout << "Error in setup_tissue: no rod cell.\n";
        return;
    }


    int number_of_rods = parameters.ints("number_of_rods");
    int N = 27;

    // Place rods in staggered rows
    // dx matches elastic confluent rest length (from built in func)
    double rad = pRodDef->phenotype.geometry.radius;
    double dx  = 2.0 * rad * 0.95238095238;     
    double row_gap = 0.86602540378 * dx;    

    for (int rod_idx = 0; rod_idx < number_of_rods; rod_idx++)
    {

        double x0 = Xmin + UniformRandom() * Xrange;
        double y0 = Ymin + UniformRandom() * Yrange;
        double z0 = 0.0;
        if (default_microenvironment_options.simulate_2D == false)
        { z0 = Zmin + UniformRandom() * Zrange; }

        // random orientation
        double theta = 2.0 * M_PI * UniformRandom();
        double ux = cos(theta), uy = sin(theta);
        double vx = -uy, vy = ux; 

        std::vector<Cell*> top(N, nullptr);
        std::vector<Cell*> bot(N, nullptr);

        for (int i = 0; i < N; i++)
        {
            double offset = (i - 0.5 * (N - 1)) * dx;

            // top row
            {
                Cell* c = create_cell(*pRodDef);
                double x = x0 + offset * ux + 0.5 * row_gap * vx;
                double y = y0 + offset * uy + 0.5 * row_gap * vy;
                c->assign_position(x, y, z0);

                // may need later 
                c->custom_data["rod_id"] = rod_idx;
                c->custom_data["row"] = 0;
                c->custom_data["col"] = i;

                top[i] = c;
            }

            // bottom row (staggered by dx/2 along axis)
            {
                Cell* c = create_cell(*pRodDef);
                double x = x0 + (offset + 0.5 * dx) * ux - 0.5 * row_gap * vx;

                double y = y0 + (offset + 0.5 * dx) * uy - 0.5 * row_gap * vy;
                c->assign_position(x, y, z0);

                c->custom_data["rod_id"] = rod_idx;
                c->custom_data["row"] = 1;
                c->custom_data["col"] = i;

                bot[i] = c;
            }
        }

        // attach along rows
        for (int i = 0; i < N - 1; i++)
        {
            attach_cells(top[i], top[i + 1]);
            attach_cells(bot[i], bot[i + 1]);
        }

        // Attach diagonally to help prevent bending 
        for (int i = 0; i < N; i++)
        {
            attach_cells(top[i], bot[i]);
            if (i > 0) attach_cells(top[i], bot[i - 1]);
        }
    }

    return;
}


void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{
    return;
}


void contact_function( Cell* pMe, Phenotype& phenoMe,
                       Cell* pOther, Phenotype& phenoOther, double dt )

{
    // only applied to cells with attachments associated with it
    standard_elastic_contact_function_confluent_rest_length(pMe, phenoMe, pOther, phenoOther,dt);
    return;
}



void tcell_update_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{
    // runs once every dt_mechanics 

    // run default mechanics 
    standard_update_cell_velocity( pCell, phenotype, dt );

    static int rod_type = get_cell_definition("rod").type;


    int k_touch = pCell->custom_data.find_variable_index("touching_rod");
    int k_time  = pCell->custom_data.find_variable_index("attached_time");
    int k_state = pCell->custom_data.find_variable_index("state");

    if( k_touch < 0 || k_time < 0 || k_state < 0 ) return; // if any of thesedon't exist 

    pCell->custom_data[k_touch] = 0.0;

    auto nearby = pCell->nearby_interacting_cells();

    for( Cell* other : nearby )
    {
        if( other->type != rod_type ) continue;

        double dx = pCell->position[0] - other->position[0];
        double dy = pCell->position[1] - other->position[1];
        double dz = pCell->position[2] - other->position[2];
        double d  = std::sqrt(dx*dx + dy*dy + dz*dz);

        double contact_dist = pCell->phenotype.geometry.radius + other->phenotype.geometry.radius;

        if( d <= contact_dist )
        {
            // T cdll is touching a rodddd
            pCell->custom_data[k_touch] = 1.0;
            break;
        }
    }

    if( pCell->custom_data[k_state] < 0.5 )
    {
        if( pCell->custom_data[k_touch] > 0.5 )
        {
            pCell->custom_data[k_time] += dt;

            const double THRESH_MIN = 30;

            if( pCell->custom_data[k_time] >= THRESH_MIN )
            {
                pCell->custom_data[k_state] = 1.0;

                #pragma omp critical
                std::cout
                << "[Tcell ACTIVATED]"
                << " id =" << pCell->ID
                << " total_time =" << pCell->custom_data[k_time]
                << " t =" << PhysiCell_globals.current_time
                << std::endl;



            }
        }
    }
}


void tcell_division_function( Cell* pParent, Cell* pDaughter )
{
    // called everytime there is a cell division, includes both parent and daughter cell 
    if( pParent == nullptr || pDaughter == nullptr ) return;

    static int t_type = get_cell_definition("T cell").type;
    if( pParent->type != t_type ) return;

    int k_state = pParent->custom_data.find_variable_index("state");
    int k_time  = pParent->custom_data.find_variable_index("attached_time");
    int k_touch = pParent->custom_data.find_variable_index("touching_rod");

    if( k_state < 0 ) return;

    // If parent is active, reset only the daughter to naive
    if( pParent->custom_data[k_state] > 0.5 )
    {
        pDaughter->custom_data[k_state] = 0.0;

        if( k_time  >= 0 ) pDaughter->custom_data[k_time]  = 0.0;
        if( k_touch >= 0 ) pDaughter->custom_data[k_touch] = 0.0;

        #pragma omp critical
        std::cerr << "[DIV] parent active -> daughter naive | parent id="
                  << pParent->ID << " daughter id=" << pDaughter->ID
                  << " t=" << PhysiCell_globals.current_time << "\n";
    }
}



std::vector<std::string> coloring_function( Cell* pCell )
{
    std::vector<std::string> output(4, "grey");

    // dead cells 
    if( pCell->phenotype.death.dead )
    {
        output[0] = "red";
        output[2] = "darkred";
        return output;
    }

    static int t_type   = get_cell_definition("T cell").type;
    static int rod_type = get_cell_definition("rod").type;

    // t cell 
    if( pCell->type == t_type )
    {
        int k_state = pCell->custom_data.find_variable_index("state");
        bool active = (k_state >= 0 && pCell->custom_data[k_state] > 0.5);

        if( active )
        {
            // active T cell
            output[0] = "#FF8C00";   // orange
            output[2] = "#FFFFFF";   // white outline
            // output[0] = "purple";
            // output[2] = "indigo";
            // output[0] = "#ff00ff";   
            // output[2] = "#aa00aa";   
        }
        else
        {
            // naive T cell
            output[0] = "blue";
            output[2] = "darkblue";
        }
    }

    // rod cells
    else if( pCell->type == rod_type )
    {
        output[0] = "green";
        output[2] = "darkgreen";
    }

    // all other cells 
    else
    {
        output[0] = "grey";
        output[2] = "darkgrey";
    }

    return output;
}

