
"""
    function update_drawio_labels_with_modifications2(drawio_file::String, csv_file::String, output_file::String;
                                               g_modifier::Function = identity,
                                               xi_modifier::Function = identity,
                                               g2_modifier::Function = identity)

Update labels in a draw.io XML file based on values from a CSV file with custom modifiers.
- Updates values associated to labels
- Formats fractions in exponential notation. Filter < 1e-15 to display as 0.0.
- When there is a property H2O = true it converts G to kgH2O/h.

```julia
update_drawio_labels_with_modifications2(
       "RUN/DEMO2018_16_CASO4/DEMO_FLOWSHEET_OFC.drawio",
       "RUN/DEMO2018_16_CASO4/prove/DEMO2018_16_CASO4-0.99975.csv",
       "RUN/DEMO2018_16_CASO4/DEMO_FLOWSHEET_OFC.drawio",
       g_modifier = x -> abs(round(x*8.314*273.15, digits=3)),
       xi_modifier = x -> x < 1e-15 ? "0.0" : @sprintf("%.2e", x),
       g2_modifier = x -> abs(round(x*(18*60*60/(3 * 1000)), digits=3)))
```
"""
function update_drawio_labels_with_modifications2(drawio_file::String, csv_file::String, output_file::String;
                                               g_modifier::Function = identity,
                                               xi_modifier::Function = identity,
                                               g2_modifier::Function = identity)
    # Read the CSV file (keep it unmodified)
    df = CSV.read(csv_file, DataFrame)
   
    # Get the component names (column names) and their values (last row)
    component_names = names(df)
    last_row_values = Vector(df[end, :])
   
    # Create a mapping from component names to their RAW values
    component_to_value = Dict(zip(component_names, last_row_values))
   
    println("Component mappings from CSV (raw values):")
    for (comp, val) in component_to_value
        println("  $comp -> $val")
    end
   
    # Parse the XML file
    doc = parse_file(drawio_file)
    root_element = LightXML.root(doc)
   
    # Find all object elements recursively
    objects = find_objects_recursive(root_element)

    # Modifier for CT
    modCT = x -> round(x, digits=2)
   
    println("\nFound $(length(objects)) objects in the diagram")
   
    updates_made = 0
   
    for obj in objects
        current_id = LightXML.attribute(obj, "id")
        current_label = LightXML.attribute(obj, "label")
        h2o_attribute = LightXML.attribute(obj, "H2O")
        
        # Handle .CT components with H2O="true"
        if endswith(current_id, ".CT") && h2o_attribute == "true"
            # Extract the base name (everything before .CT)
            base_name = current_id[1:end-3]  # Remove ".CT"
            
            # Look for corresponding .G and .xT components
            g_component = base_name * ".G"
            xt_component = base_name * ".xT"
            
            if haskey(component_to_value, g_component) && haskey(component_to_value, xt_component)
                # Get the raw values
                raw_g_value = component_to_value[g_component]
                raw_xt_value = component_to_value[xt_component]
                
                # Convert to numbers if needed
                g_value = isa(raw_g_value, String) ? parse(Float64, raw_g_value) : Float64(raw_g_value)
                xt_value = isa(raw_xt_value, String) ? parse(Float64, raw_xt_value) : Float64(raw_xt_value)
                
                # Calculate CT value (hardcoded calculation)
                new_value = modCT(g_value * xt_value * 3 * 9619.0 / (g_value*18/(3*1000)))  # Replace this line with your specific calculation
                
                println("Updating '$current_id' (.CT component with H2O=true): '$current_label' -> '$new_value' (calculated from G=$g_value, xT=$xt_value)")
                
                # Update the label
                LightXML.set_attribute(obj, "label", string(new_value))
                updates_made += 1
            else
                println("Warning: Could not find corresponding .G or .xT component for '$current_id'")
            end
            
        # Check if the ID matches any component name (existing logic)
        elseif haskey(component_to_value, current_id)
            raw_value = component_to_value[current_id]
           
            # Convert string to number if needed
            num_value = isa(raw_value, String) ? parse(Float64, raw_value) : Float64(raw_value)
           
            # Apply the appropriate modifier based on component type
            if endswith(current_id, ".G")
                # For .G components, check if H2O="true"
                if h2o_attribute == "true"
                    new_value = g2_modifier(num_value)
                    println("Updating '$current_id' (.G component with H2O=true): '$current_label' -> '$new_value'")
                else
                    new_value = g_modifier(num_value)
                    println("Updating '$current_id' (.G component): '$current_label' -> '$new_value'")
                end
            elseif occursin(r"\.x.+$", current_id)
                # Apply xi_modifier for .xi components
                new_value = xi_modifier(num_value)
                println("Updating '$current_id' (.xi component): '$current_label' -> '$new_value'")
            else
                # Use raw value for other components
                new_value = raw_value
                println("Updating '$current_id': '$current_label' -> '$new_value'")
            end
           
            # Update the label
            LightXML.set_attribute(obj, "label", string(new_value))
            updates_made += 1
        end
    end
   
    println("\nMade $updates_made updates")
   
    # Save the modified XML
    LightXML.save_file(doc, output_file)
    LightXML.free(doc)
   
    println("Updated diagram saved to: $output_file")
end


"""
    set_drawio_labels_to_na(old_drawio_file, new_drawio_file)

Function to set all albels to "NA"

## Example usage:
```julia
set_drawio_labels_to_na("RUN/DEMO2018_16_CASO4/DEMO_FLOWSHEET_OFC.drawio",
                        "RUN/DEMO2018_16_CASO4/DEMO_FLOWSHEET_OFC.drawio")
```
"""
function set_drawio_labels_to_na(drawio_file::String, output_file::String)
    # Parse the XML file
    doc = parse_file(drawio_file)
    root_element = LightXML.root(doc)
    
    # Find all object elements recursively
    objects = find_objects_recursive(root_element)
    
    println("Found $(length(objects)) objects in the diagram")
    
    updates_made = 0
    
    for obj in objects
        current_id = LightXML.attribute(obj, "id")
        current_label = LightXML.attribute(obj, "label")
        
        # Check if ID matches the patterns we want to set to "NA"
        should_set_na = false
        
        if current_id !== nothing
            # Check if ID ends with ".G"
            if endswith(current_id, ".G")
                should_set_na = true
                println("Setting '$current_id' to NA (ends with .G): '$current_label' -> 'NA'")
            # Check if ID matches pattern ".x" + any character (e.g., ".x1", ".x2", ".xa", etc.)
            elseif occursin(r"\.x.+$", current_id)
                should_set_na = true
                println("Setting '$current_id' to NA (matches .x pattern): '$current_label' -> 'NA'")
            end
        end
        
        if should_set_na
            # Update the label to "NA"
            LightXML.set_attribute(obj, "label", "NA")
            updates_made += 1
        end
    end
    
    println("\nMade $updates_made updates to 'NA'")
    
    # Save the modified XML
    LightXML.save_file(doc, output_file)
    LightXML.free(doc)
    
    println("Updated diagram saved to: $output_file")
end


"""
    modify_csv_values(df::DataFrame, g_modifier::Function, xi_modifier::Function)
Modifies CSV values before updating XML
"""
function modify_csv_values(df::DataFrame, g_modifier::Function, xi_modifier::Function)
    # Make a copy to avoid modifying the original
    modified_df = copy(df)
    
    column_names = names(modified_df)
    
    println("Modifying CSV values:")
    
    for col_name in column_names
        if endswith(col_name, ".G")
            # Modify .G columns
            println("  Modifying column '$col_name' (ends with .G)")
            modified_df[!, col_name] = g_modifier.(modified_df[!, col_name])
        elseif occursin(r"\.x.+$", col_name)  # Regex to match .x followed by any character
            # Modify .xi columns (where i is any character)
            println("  Modifying column '$col_name' (ends with .x + any character)")
            modified_df[!, col_name] = xi_modifier.(modified_df[!, col_name])
        end
    end
    
    return modified_df
end



"""
    find_objects_recursive(element)
Helper function to find all object elements in a drawio diagram recursively
"""
function find_objects_recursive(element)
    objects = []
    
    # Check if current element is an object
    if LightXML.name(element) == "object"
        push!(objects, element)
    end
    
    # Recursively search children
    for child in LightXML.child_elements(element)
        append!(objects, find_objects_recursive(child))
    end
    
    return objects
end