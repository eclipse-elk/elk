/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.data;

import java.util.LinkedList;
import java.util.List;

/**
 * Data type used to store information for a layout category. A layout category defines the basic
 * graph drawing approach employed by a layout algorithm, for instance the layer-based approach
 * or the force-based approach. The category of an algorithm can be left empty, in which case it
 * is displayed as "Other" in the user interface.
 *
 * @author msp
 */
public final class LayoutCategoryData implements ILayoutMetaData {
    
    /** identifier of the layout type. */
    private final String id;
    /** user friendly name of the layout type. */
    private final String name;
    /** detail description. */
    private final String description;
    /** the list of layout algorithms that are registered for this category. */
    private final List<LayoutAlgorithmData> layouters = new LinkedList<LayoutAlgorithmData>();
    
    /**
     * Create a layout category data entry.
     */
    private LayoutCategoryData(final Builder builder) {
        this.id = builder.id;
        this.name = builder.name;
        this.description = builder.description;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof LayoutCategoryData) {
            return this.id.equals(((LayoutCategoryData) obj).id);
        }
        return false;
    }
    
    @Override
    public int hashCode() {
        return id.hashCode();
    }
    
    @Override
    public String toString() {
        return "Layout Type: " + id;
    }
    
    /**
     * Returns the list of layout algorithms that are registered for this type.
     * 
     * @return the layouters
     */
    public List<LayoutAlgorithmData> getLayouters() {
        return layouters;
    }

    /**
     * Returns the layout type identifier.
     * 
     * @return the layout type identifier
     */
    public String getId() {
        return id;
    }
    
    /**
     * Returns the name of the layout type.
     * 
     * @return the name
     */
    public String getName() {
        return name;
    }
    
    /**
     * Returns the description.
     * 
     * @return the description
     */
    public String getDescription() {
        return description;
    }
    
    /**
     * Builder for {@link LayoutCategoryData} instances.
     */
    public static class Builder {
        
        private String id;
        private String name;
        private String description;
        
        /**
         * Create an instance with the configured values.
         */
        public LayoutCategoryData create() {
            return new LayoutCategoryData(this);
        }
        
        /**
         * Configure the {@link LayoutCategoryData#getId() id}.
         */
        public Builder id(final String aid) {
            this.id = aid;
            return this;
        }
        
        /**
         * Configure the {@link LayoutCategoryData#getName() name}.
         */
        public Builder name(final String aname) {
            this.name = aname;
            return this;
        }
        
        /**
         * Configure the {@link LayoutCategoryData#getDescription() description}.
         */
        public Builder description(final String adescription) {
            this.description = adescription;
            return this;
        }
        
    }
    
}
