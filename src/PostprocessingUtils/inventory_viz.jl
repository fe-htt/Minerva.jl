

"""
    update_drawio_inventory_colours(drawio_file::String, csv_file::String, output_file::String;
                               inventory_column::String = "NT",
                               colour_scheme::Symbol = :viridis,
                               log_scale::Bool = true,
                               Vplasma::Union{Float64, Nothing} = nothing,
                               linear_threshold::Float64 = 1.0,
                               nt_scaling_factor::Float64 = 1.0)

Updates draw.io diagram component colours based on inventory data from CSV.

# Example usage
```julia
update_drawio_inventory_colours(
    "RUN/DEMO2018_16_CASO4/DEMO_inventories_custom.drawio", 
    "RUN/DEMO2018_16_CASO4/prove/DEMO2018_16_CASO4-0.99975.csv",
    "RUN/DEMO2018_16_CASO4/DEMO_inventories_custom.drawio",
    inventory_column = "NT",
    linear_threshold = 1.0,  # Values above 1.0 mol use linear scaling
    #log_scale = false,
    colour_scheme = :jet,  # Options: :viridis, :plasma, :inferno, :cividis, :thermal, :simple, :jet, :turbo
    Vplasma = 2579.0,
    nt_scaling_factor = 3.0  # scale all NT values by this factor
)
```
"""
function update_drawio_inventory_colours(drawio_file::String, csv_file::String, output_file::String;
                                       inventory_column::String = "NT",
                                       colour_scheme::Symbol = :viridis,
                                       log_scale::Bool = true,  # Keep as default but now means "hybrid" scaling
                                       Vplasma::Union{Float64, Nothing} = nothing,
                                       linear_threshold::Float64 = 1.0,  # Values above this use linear scaling
                                       nt_scaling_factor::Float64 = 1.0)  # NEW: scaling factor for NT values
    
    # Read the CSV file
    df = CSV.read(csv_file, DataFrame)
    
    # Get component names (column names) and their values (last row)
    component_names = names(df)
    last_row_values = Vector(df[end, :])
    
    # Create mapping from component names to inventory values
    component_to_inventory = Dict(zip(component_names, last_row_values))
    
    # Filter for components that have inventory data
    inventory_components = filter(p -> occursin(inventory_column, first(p)), component_to_inventory)
    
    if isempty(inventory_components)
        error("No components found with inventory column '$inventory_column'")
    end
    
    # Extract base component names (remove .NT suffix) and handle special cases
    base_components = Dict{String, Float64}()
    for (comp_name, value) in inventory_components
        base_name = replace(comp_name, ".$(inventory_column)" => "")
        numeric_value = isa(value, String) ? parse(Float64, value) : Float64(value)
        
        # Apply scaling factor to NT values
        scaled_value = numeric_value * nt_scaling_factor
        base_components[base_name] = scaled_value
        
        # Log the scaling application if factor is not 1.0
        if nt_scaling_factor != 1.0
            println("Scaled $(comp_name): $(Printf.@sprintf("%.3e", numeric_value)) → $(Printf.@sprintf("%.3e", scaled_value)) mol (factor: $(nt_scaling_factor))")
        end
    end
    
    # Handle special case for myTorus: calculate NT as nT * Vplasma  
    if Vplasma !== nothing && haskey(component_to_inventory, "myTorus.nT")
        nT_value = component_to_inventory["myTorus.nT"]
        
        # Convert to numeric if needed
        nT_numeric = isa(nT_value, String) ? parse(Float64, nT_value) : Float64(nT_value)
        
        # Calculate tritium inventory for torus
        torus_NT = nT_numeric * Vplasma
        
        # Apply scaling factor to torus inventory
        scaled_torus_NT = torus_NT * nt_scaling_factor
        base_components["myTorus"] = scaled_torus_NT
        
        if nt_scaling_factor != 1.0
            println("Special calculation for myTorus with scaling:")
            println("  Original: NT = nT × Vplasma = $(Printf.@sprintf("%.3e", nT_numeric)) × $(Printf.@sprintf("%.3e", Vplasma)) = $(Printf.@sprintf("%.3e", torus_NT)) mol")
            println("  Scaled: NT = $(Printf.@sprintf("%.3e", torus_NT)) × $(nt_scaling_factor) = $(Printf.@sprintf("%.3e", scaled_torus_NT)) mol")
        else
            println("Special calculation for myTorus: NT = nT × Vplasma = $(Printf.@sprintf("%.3e", nT_numeric)) × $(Printf.@sprintf("%.3e", Vplasma)) = $(Printf.@sprintf("%.3e", scaled_torus_NT)) mol")
        end
    elseif Vplasma === nothing && haskey(component_to_inventory, "myTorus.nT")
        println("Warning: Found myTorus.nT but Vplasma parameter not provided, skipping torus inventory calculation")
    elseif Vplasma !== nothing && !haskey(component_to_inventory, "myTorus.nT")
        println("Warning: Vplasma parameter provided but myTorus.nT not found in data, skipping torus inventory calculation")
    end
    
    # Display scaling factor information
    if nt_scaling_factor != 1.0
        println("\n" * "="^60)
        println("NT SCALING APPLIED")
        println("="^60)
        println("Scaling factor: $(nt_scaling_factor)")
        println("All NT inventory values have been multiplied by this factor")
        println("This affects visualization only - original data unchanged")
        println("="^60)
    end
    
    # Parse the XML file first to get available diagram components
    doc = parse_file(drawio_file)
    root_element = LightXML.root(doc)
    
    # Find all object elements recursively
    objects = find_objects_recursive(root_element)
    println("\nFound $(length(objects)) objects in the diagram")
    
    # Get list of component IDs that are actually in the diagram
    diagram_component_ids = Set{String}()
    for obj in objects
        current_id = LightXML.attribute(obj, "id")
        if current_id !== nothing
            push!(diagram_component_ids, current_id)
        end
    end
    
    println("Diagram component IDs: $(sort(collect(diagram_component_ids)))")
    
    # Filter base_components to only include those that exist in the diagram
    diagram_components = Dict{String, Float64}()
    for (comp_name, inventory) in base_components
        if comp_name in diagram_component_ids
            diagram_components[comp_name] = inventory
        end
    end
    
    if isempty(diagram_components)
        println("Warning: No components with inventory data found in the diagram")
        println("Available inventory components: $(sort(collect(keys(base_components))))")
        println("Diagram component IDs: $(sort(collect(diagram_component_ids)))")
        return
    end
    
    println("Found $(length(diagram_components)) components with both inventory data and diagram representation:")
    for (comp, inv) in sort(collect(diagram_components), by=x->x[2], rev=true)
        println("  $comp: $(Printf.@sprintf("%.3e", inv)) mol")
    end
    
    # Calculate inventory range for colour scaling (only for diagram components)
    inventories = collect(values(diagram_components))
    
    if log_scale
        # Hybrid log-linear scaling approach
        positive_inventories = filter(x -> x > 0, inventories)
        if isempty(positive_inventories)
            error("No positive inventory values found for hybrid scaling")
        end
        
        min_inv = minimum(positive_inventories)
        max_inv = maximum(positive_inventories)
        
        # Separate inventories into low (log-scaled) and high (linear-scaled) ranges
        low_inventories = filter(x -> x <= linear_threshold, positive_inventories)
        high_inventories = filter(x -> x > linear_threshold, positive_inventories)
        
        println("\nHybrid scaling analysis:")
        println("  Total range: $(Printf.@sprintf("%.3e", min_inv)) to $(Printf.@sprintf("%.3e", max_inv)) mol")
        if nt_scaling_factor != 1.0
            println("  (After applying NT scaling factor: $(nt_scaling_factor))")
        end
        println("  Linear threshold: $(linear_threshold) mol")
        println("  Low values (log-scaled): $(length(low_inventories)) components")
        println("  High values (linear-scaled): $(length(high_inventories)) components")
        
        if !isempty(low_inventories)
            println("  Low range: $(Printf.@sprintf("%.3e", minimum(low_inventories))) to $(Printf.@sprintf("%.3e", maximum(low_inventories))) mol")
        end
        if !isempty(high_inventories)
            println("  High range: $(Printf.@sprintf("%.3e", minimum(high_inventories))) to $(Printf.@sprintf("%.3e", max_inv)) mol")
        end
        
        # Calculate percentiles for the high-value range
        if !isempty(high_inventories)
            sorted_high = sort(high_inventories)
            p25 = sorted_high[max(1, round(Int, 0.25 * length(sorted_high)))]
            p50 = sorted_high[max(1, round(Int, 0.50 * length(sorted_high)))]
            p75 = sorted_high[max(1, round(Int, 0.75 * length(sorted_high)))]
            println("  High range percentiles:")
            println("    25th: $(Printf.@sprintf("%.2f", p25)) mol")
            println("    50th: $(Printf.@sprintf("%.2f", p50)) mol") 
            println("    75th: $(Printf.@sprintf("%.2f", p75)) mol")
        end
        
    else
        # Pure linear scaling
        min_inv = minimum(inventories)
        max_inv = maximum(inventories)
        range_info = "Inventory range (linear): $(Printf.@sprintf("%.3e", min_inv)) to $(Printf.@sprintf("%.3e", max_inv)) mol"
        if nt_scaling_factor != 1.0
            range_info *= " (After NT scaling factor: $(nt_scaling_factor))"
        end
        println("\n" * range_info)
    end
    
    # Handle case where all inventories are the same
    if min_inv == max_inv
        println("Warning: All inventory values are identical ($(Printf.@sprintf("%.3e", min_inv)) mol)")
        println("All components will receive the same colour")
        inv_range = 1.0  # Avoid division by zero
    else
        inv_range = max_inv - min_inv
    end
    
    updates_made = 0
    
    for obj in objects
        current_id = LightXML.attribute(obj, "id")
        
        if current_id !== nothing && haskey(diagram_components, current_id)
            inventory = diagram_components[current_id]
            
            # Calculate normalised value using hybrid log-linear scaling
            if min_inv == max_inv
                norm_value = 0.5  # Middle colour when all values are the same
            elseif log_scale && inventory > 0
                # Hybrid approach: log for low values, linear for high values
                if inventory <= linear_threshold
                    # Log scaling for low values (good for 10^-4 to 1.0 range)
                    low_min = min_inv
                    low_max = min(linear_threshold, max_inv)
                    
                    if low_min == low_max
                        norm_value = 0.0  # All low values get minimum colour
                    else
                        log_norm = (log10(inventory) - log10(low_min)) / (log10(low_max) - log10(low_min))
                        norm_value = log_norm * 0.3  # Use bottom 30% of colour range for low values
                    end
                else
                    # Linear scaling for high values (good for 1.0 to 100+ range)  
                    high_min = linear_threshold
                    high_max = max_inv
                    
                    if high_min == high_max
                        norm_value = 0.7  # Default position if no range
                    else
                        linear_norm = (inventory - high_min) / (high_max - high_min)
                        norm_value = 0.3 + linear_norm * 0.7  # Use top 70% of colour range for high values
                    end
                end
            else
                # Pure linear scaling
                norm_value = (inventory - min_inv) / inv_range
            end
            
            # Clamp to [0, 1] range
            norm_value = clamp(norm_value, 0.0, 1.0)
            
            # Generate colour based on scheme
            colour = get_colour_from_scheme(norm_value, colour_scheme)
            
            # Convert to hex format for draw.io (no alpha for better visibility)
            hex_colour = string("#", hex(colour))
            
            # MODIFIED: Format inventory value for label with 0.01 threshold for exponential notation
            if inventory >= 1e3 || inventory < 0.01
                inventory_label = @sprintf("%.1e", inventory)
            else
                inventory_label = @sprintf("%.2f", inventory)
            end
            
            # Update object label attribute
            LightXML.set_attribute(obj, "label", inventory_label)
            
            # Find the mxCell child element and update its style
            mxcell = nothing
            for child in LightXML.child_elements(obj)
                if LightXML.name(child) == "mxCell"
                    mxcell = child
                    break
                end
            end
            
            if mxcell !== nothing
                # Update fillColor in mxCell style attribute (remove alpha for better visibility)
                current_style = LightXML.attribute(mxcell, "style")
                if current_style !== nothing
                    # Replace existing fillColor
                    new_style = update_fill_colour(current_style, hex_colour)
                else
                    # Create new style with fillColor
                    new_style = "fillColor=$(hex_colour);strokeColor=#000000;"
                end
                
                LightXML.set_attribute(mxcell, "style", new_style)
            else
                println("Warning: No mxCell found for object '$current_id'")
            end
            
            println("Updated '$current_id': inventory=$(inventory_label) mol, colour=$hex_colour")
            updates_made += 1
        end
    end
    
    println("\nMade $updates_made colour updates")
    
    # Update legend if present
    update_legend(objects, min_inv, max_inv, colour_scheme, log_scale, inv_range, linear_threshold)
    
    # Save the modified XML
    LightXML.save_file(doc, output_file)
    LightXML.free(doc)
    
    println("Coloured diagram saved to: $output_file")
    
    # Generate legend information
    generate_colour_legend(min_inv, max_inv, colour_scheme, log_scale, nt_scaling_factor)
end


"""
    get_colour_from_scheme(norm_value::Float64, scheme::Symbol)
Generate colour from normalised value using specified colour scheme

## Available schemes:
- :viridis
- :plasma
- :inferno
- :jet
- :turbo
- :thermal
- :simple
default is :viridis
"""
function get_colour_from_scheme(norm_value::Float64, scheme::Symbol)
    
    # Remove additional power scaling since we now use sigmoid in the main function
    # The sigmoid already provides the adaptive distribution we need
    
    if scheme == :viridis
        return Colors.RGB(viridis_colour(norm_value)...)
    elseif scheme == :plasma
        return Colors.RGB(plasma_colour(norm_value)...)
    elseif scheme == :inferno
        return Colors.RGB(inferno_colour(norm_value)...)
    elseif scheme == :jet
        return Colors.RGB(jet_colour(norm_value)...)
    elseif scheme == :turbo
        return Colors.RGB(turbo_colour(norm_value)...)
    elseif scheme == :thermal
        # Simple blue to red thermal scale
        return Colors.RGB(norm_value, 0.0, 1.0 - norm_value)
    elseif scheme == :simple
        # High contrast simple scheme for clear distinction
        if norm_value < 0.2
            return Colors.RGB(0.0, 0.0, 1.0)  # Blue for lowest
        elseif norm_value < 0.4
            return Colors.RGB(0.0, 0.5, 1.0)  # Light blue
        elseif norm_value < 0.6
            return Colors.RGB(0.0, 1.0, 0.5)  # Green
        elseif norm_value < 0.8
            return Colors.RGB(1.0, 1.0, 0.0)  # Yellow
        else
            return Colors.RGB(1.0, 0.0, 0.0)  # Red for highest
        end
    else
        # Default viridis
        return Colors.RGB(viridis_colour(norm_value)...)
    end
end


"""
    jet_colour(t::Float64)
Jet colourmap implementation - classic vibrant blue-cyan-green-yellow-red
"""
function jet_colour(t::Float64)
    t = clamp(t, 0.0, 1.0)
    
    # Classic jet colourmap with proper vibrant colours
    if t < 0.25
        # Blue to Cyan
        r = 0.0
        g = 4.0 * t
        b = 1.0
    elseif t < 0.5
        # Cyan to Green
        r = 0.0
        g = 1.0
        b = 1.0 - 4.0 * (t - 0.25)
    elseif t < 0.75
        # Green to Yellow
        r = 4.0 * (t - 0.5)
        g = 1.0
        b = 0.0
    else
        # Yellow to Red
        r = 1.0
        g = 1.0 - 4.0 * (t - 0.75)
        b = 0.0
    end
    
    return (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
end


"""
    turbo_colour(t::Float64)
Turbo colourmap implementation
"""
function turbo_colour(t::Float64)
    t = clamp(t, 0.0, 1.0)
    
    # Turbo colourmap - polynomial approximation for vibrant colours
    # Red channel
    r = 0.18995 + t * (4.55940 + t * (-42.3277 + t * (130.5 + t * (-118.1 + t * 35.0))))
    
    # Green channel  
    g = 0.07176 + t * (2.81160 + t * (-2.85120 + t * (-18.2400 + t * (52.1600 + t * (-40.0)))))
    
    # Blue channel
    b = 0.23217 + t * (8.24020 + t * (-60.1920 + t * (109.600 + t * (-88.7500 + t * 25.0))))
    
    return (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
end


"""
    viridis_colour(t::Float64)
Viridis colourmap implementation
"""
function viridis_colour(t::Float64)
    t = clamp(t, 0.0, 1.0)
    c0 = [0.2777273272234177, 0.005407344544966578, 0.3340998053353061]
    c1 = [0.1050930431085774, 1.404613529898575, 1.384590162594685]
    c2 = [-0.3308618287255563, 0.214847559468213, 0.09509516302823659]
    c3 = [-4.634230498983486, -5.799100973351585, -19.33244095627987]
    c4 = [6.228269936347081, 14.17993336680509, 56.69055260068105]
    c5 = [4.776384997670288, -13.74514537774601, -65.35303263337234]
    c6 = [-5.435455855934631, 4.645852612178535, 26.3124352495832]
    
    r = c0[1] + c1[1]*t + c2[1]*t^2 + c3[1]*t^3 + c4[1]*t^4 + c5[1]*t^5 + c6[1]*t^6
    g = c0[2] + c1[2]*t + c2[2]*t^2 + c3[2]*t^3 + c4[2]*t^4 + c5[2]*t^5 + c6[2]*t^6
    b = c0[3] + c1[3]*t + c2[3]*t^2 + c3[3]*t^3 + c4[3]*t^4 + c5[3]*t^5 + c6[3]*t^6
    
    return (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
end


"""
    plasma_colour(t::Float64)
Plasma colourmap implementation
"""
function plasma_colour(t::Float64)
    t = clamp(t, 0.0, 1.0)
    c0 = [0.05873234392399702, 0.02333670892565664, 0.5433401826748754]
    c1 = [2.176514634195958, 4.478572498805903, 2.053257277540509]
    c2 = [-2.689460476458034, -3.5498806451913123, -2.1974017513702324]
    c3 = [6.130348345893603, 5.3639202545340515, -7.1455098160204425]
    c4 = [11.11386993844262, -3.2455678110211893, 24.549673922988808]
    c5 = [-41.786737506084216, 22.93153465461149, -56.16207193982797]
    c6 = [77.43463343603947, -60.51342071671263, 67.52394027048522]
    c7 = [-71.3175974177072, 67.4984722895303, -36.743181732963204]
    
    r = c0[1] + c1[1]*t + c2[1]*t^2 + c3[1]*t^3 + c4[1]*t^4 + c5[1]*t^5 + c6[1]*t^6 + c7[1]*t^7
    g = c0[2] + c1[2]*t + c2[2]*t^2 + c3[2]*t^3 + c4[2]*t^4 + c5[2]*t^5 + c6[2]*t^6 + c7[2]*t^7
    b = c0[3] + c1[3]*t + c2[3]*t^2 + c3[3]*t^3 + c4[3]*t^4 + c5[3]*t^5 + c6[3]*t^6 + c7[3]*t^7
    
    return (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
end


"""
    inferno_colour(t::Float64)
Inferno colourmap implementation
"""
function inferno_colour(t::Float64)
    t = clamp(t, 0.0, 1.0)
    c0 = [0.0002189403691192265, 0.001651004631001012, -0.01948089843709184]
    c1 = [0.1065134194856116, 0.5639564367884091, 3.932712388889277]
    c2 = [11.60249308247187, -3.972853965665698, -15.9423941062914]
    c3 = [-41.70399613139459, 17.43639888205313, 44.35414519872813]
    c4 = [77.162935699427, -33.40235894210092, -81.80730925738993]
    c5 = [-71.31942824499214, 32.62606426397723, 73.20951985803202]
    c6 = [25.13112622477341, -12.24266895238567, -23.07032500287172]
    
    r = c0[1] + c1[1]*t + c2[1]*t^2 + c3[1]*t^3 + c4[1]*t^4 + c5[1]*t^5 + c6[1]*t^6
    g = c0[2] + c1[2]*t + c2[2]*t^2 + c3[2]*t^3 + c4[2]*t^4 + c5[2]*t^5 + c6[2]*t^6
    b = c0[3] + c1[3]*t + c2[3]*t^2 + c3[3]*t^3 + c4[3]*t^4 + c5[3]*t^5 + c6[3]*t^6
    
    return (clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0))
end


"""
    update_fill_colour(style::String, hex_colour::String)
Update or add fillColor in draw.io style string
"""
function update_fill_colour(style::String, hex_colour::String)
    # Remove existing fillColor if present
    style_clean = replace(style, r"fillColor=[^;]*;" => "")
    
    # Add new fillColor at the beginning
    return "fillColor=$(hex_colour);" * style_clean
end


"""
    generate_colour_legend(min_inv::Float64, max_inv::Float64, scheme::Symbol, log_scale::Bool, nt_scaling_factor::Float64 = 1.0)
Generate legend information for the colour scheme
"""
function generate_colour_legend(min_inv::Float64, max_inv::Float64, scheme::Symbol, log_scale::Bool, nt_scaling_factor::Float64 = 1.0)
    println("\n" * "="^60)
    println("COLOUR LEGEND")
    println("="^60)
    println("Colour scheme: $scheme")
    println("Scale type: $(log_scale ? "logarithmic" : "linear")")
    println("Range: $(Printf.@sprintf("%.3e", min_inv)) to $(Printf.@sprintf("%.3e", max_inv)) mol")
    if nt_scaling_factor != 1.0
        println("NT scaling factor applied: $(nt_scaling_factor)")
    end
    println("\nColour mapping (normalised values):")
    
    for norm_val in [0.0, 0.25, 0.5, 0.75, 1.0]
        if log_scale
            actual_val = 10^(log10(min_inv) + norm_val * (log10(max_inv) - log10(min_inv)))
        else
            actual_val = min_inv + norm_val * (max_inv - min_inv)
        end
        colour = get_colour_from_scheme(norm_val, scheme)
        hex_colour = string("#", hex(colour))
        println("  $(Printf.@sprintf("%3.0f", norm_val*100))%: $(Printf.@sprintf("%.3e", actual_val)) mol → $hex_colour")
    end
    println("="^60)
end


"""
    update_legend(objects, min_inv::Float64, max_inv::Float64, colour_scheme::Symbol, log_scale::Bool, inv_range::Float64, linear_threshold::Float64 = 1.0)
Update segmented legend gradient and labels to match hybrid scaling used for components
"""
function update_legend(objects, min_inv::Float64, max_inv::Float64, colour_scheme::Symbol, log_scale::Bool, inv_range::Float64, linear_threshold::Float64 = 1.0)
    # Find legend objects - only looking for segments and labels
    legend_segments = Dict{String, Any}()
    legend_labels = Dict{String, Any}()
    
    # Debug: Print all object IDs that start with "legend"
    println("Debug: All objects with IDs starting with 'legend':")
    for obj in objects
        current_id = LightXML.attribute(obj, "id")
        if current_id !== nothing && startswith(current_id, "legend")
            println("  Found object ID: '$current_id'")
            
            if occursin("-", current_id)
                # This is a legend segment (e.g., "legend_0-25")
                legend_segments[current_id] = obj
                println("    → Classified as segment")
            else
                # This is a legend label (e.g., "legend_0")
                legend_labels[current_id] = obj
                println("    → Classified as label")
            end
        end
    end
    
    println("Found legend segments: $(sort(collect(keys(legend_segments))))")
    println("Found legend labels: $(sort(collect(keys(legend_labels))))")
    
    # Define the 5 key points for color continuity (0%, 25%, 50%, 75%, 100%)
    legend_points = [0.0, 0.25, 0.5, 0.75, 1.0]
    segment_colors = []
    
    # Calculate colors at each point using the same hybrid scaling as components
    for norm_val in legend_points
        # Apply the same hybrid scaling logic as components
        if log_scale && min_inv != max_inv  # Need to check if we have valid inventory range
            if norm_val <= 0.3
                # Low range (log-scaled region)
                low_min = min_inv
                low_max = min(linear_threshold, max_inv)
                if low_min == low_max
                    component_norm = 0.0
                else
                    log_interp = norm_val / 0.3
                    # Calculate the actual component norm value that would be used
                    # This simulates the same calculation as in the main component loop
                    simulated_inventory = 10^(log10(low_min) + log_interp * (log10(low_max) - log10(low_min)))
                    # Apply same hybrid scaling as components
                    log_norm_sim = (log10(simulated_inventory) - log10(low_min)) / (log10(low_max) - log10(low_min))
                    component_norm = log_norm_sim * 0.3
                end
            else
                # High range (linear-scaled region)
                high_min = linear_threshold
                high_max = max_inv
                if high_min == high_max
                    component_norm = 0.7
                else
                    linear_interp = (norm_val - 0.3) / 0.7
                    # Calculate the actual component norm value that would be used
                    # This simulates the same calculation as in the main component loop
                    simulated_inventory = high_min + linear_interp * (high_max - high_min)
                    linear_norm_sim = (simulated_inventory - high_min) / (high_max - high_min)
                    component_norm = 0.3 + linear_norm_sim * 0.7
                end
            end
        else
            component_norm = norm_val
        end
        
        # Get the actual color that would be used for components
        color = get_colour_from_scheme(component_norm, colour_scheme)
        push!(segment_colors, color)
        println("Legend point $(norm_val*100)%: component_norm=$(component_norm), color=#$(hex(color))")
    end
    
    # Update each legend segment with proper gradient
    segment_names = ["legend_0-25", "legend_25-50", "legend_50-75", "legend_75-100"]
    
    for (i, segment_name) in enumerate(segment_names)
        println("\n=== Processing segment: '$segment_name' ===")
        if haskey(legend_segments, segment_name)
            segment_obj = legend_segments[segment_name]
            println("✓ Found segment object")
            
            # Get start and end colors for this segment
            start_color = segment_colors[i]
            end_color = segment_colors[i + 1]
            
            start_hex = string("#", hex(start_color))
            end_hex = string("#", hex(end_color))
            
            println("✓ Segment colors: $(start_hex) → $(end_hex)")
            
            # Find the mxCell and update its gradient
            mxcell_found = false
            children_count = 0
            
            for child in LightXML.child_elements(segment_obj)
                children_count += 1
                child_name = LightXML.name(child)
                println("  Child $(children_count): <$(child_name)>")
                
                if child_name == "mxCell"
                    mxcell_found = true
                    println("  ✓ Found mxCell child")
                    
                    current_style = LightXML.attribute(child, "style")
                    println("  Current style: '$(current_style)'")
                    
                    # Create new style with gradient
                    if current_style !== nothing
                        # Remove existing fillColor and gradientColor
                        new_style = replace(current_style, r"fillColor=[^;]*;" => "")
                        new_style = replace(new_style, r"gradientColor=[^;]*;" => "")
                        
                        # Add new colours at the beginning
                        new_style = "fillColor=$(start_hex);gradientColor=$(end_hex);" * new_style
                        
                        # Ensure gradientDirection is set
                        if !occursin("gradientDirection", new_style)
                            new_style = new_style * "gradientDirection=east;"
                        end
                    else
                        new_style = "fillColor=$(start_hex);gradientColor=$(end_hex);gradientDirection=east;"
                    end
                    
                    println("  New style: '$(new_style)'")
                    
                    # Actually set the attribute
                    try
                        LightXML.set_attribute(child, "style", new_style)
                        println("  ✓ Successfully updated style attribute")
                        
                        # Verify the update
                        verification_style = LightXML.attribute(child, "style")
                        println("  Verification - style now: '$(verification_style)'")
                        
                    catch e
                        println("  ✗ ERROR setting style attribute: $(e)")
                    end
                    
                    break
                end
            end
            
            println("  Total children found: $(children_count)")
            if !mxcell_found
                println("  ✗ WARNING: No mxCell found in segment $(segment_name)")
                println("  Available children: $(join([LightXML.name(child) for child in LightXML.child_elements(segment_obj)], ", "))")
            else
                println("  ✓ mxCell found and processed")
            end
        else
            println("✗ Legend segment '$(segment_name)' not found in collected segments")
            println("  Available segments: $(sort(collect(keys(legend_segments))))")
        end
        println("=== End segment $(segment_name) ===\n")
    end
    
    # Update legend labels to match hybrid scaling with MODIFIED formatting for 0.01 threshold
    legend_percentages = ["0", "25", "50", "75", "100"]
    legend_values = [0.0, 0.25, 0.5, 0.75, 1.0]
    
    for (i, percentage) in enumerate(legend_percentages)
        label_id = "legend_$(percentage)"
        if haskey(legend_labels, label_id)
            norm_val = legend_values[i]
            
            # Calculate actual inventory value using SAME hybrid scaling as components
            if min_inv == max_inv
                actual_val = min_inv
            elseif log_scale
                # Use the same hybrid logic as in the main coloring function
                if norm_val <= 0.3
                    # Low range (log-scaled region)
                    low_min = min_inv
                    low_max = min(linear_threshold, max_inv)
                    if low_min == low_max
                        actual_val = low_min
                    else
                        log_interp = norm_val / 0.3  # Rescale to [0,1] for low range
                        log_val = log10(low_min) + log_interp * (log10(low_max) - log10(low_min))
                        actual_val = 10^log_val
                    end
                else
                    # High range (linear-scaled region)
                    high_min = linear_threshold
                    high_max = max_inv
                    if high_min == high_max
                        actual_val = high_min
                    else
                        linear_interp = (norm_val - 0.3) / 0.7  # Rescale to [0,1] for high range
                        actual_val = high_min + linear_interp * (high_max - high_min)
                    end
                end
            else
                actual_val = min_inv + norm_val * inv_range
            end
            
            # MODIFIED: Format the label consistently with new 0.01 threshold
            if actual_val >= 1e3 || actual_val < 0.01
                label_text = @sprintf("%.1e", actual_val)
            else
                label_text = @sprintf("%.2f", actual_val)
            end
            
            # Update the label
            LightXML.set_attribute(legend_labels[label_id], "label", label_text)
            println("Updated $(label_id): $(label_text) mol")
        else
            println("Warning: Legend label '$(label_id)' not found")
        end
    end
    
    println("Segmented legend update completed")
end


