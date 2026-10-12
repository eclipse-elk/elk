/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
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

import java.util.EnumSet;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.util.IFactory;
import org.eclipse.elk.core.util.InstancePool;
import org.eclipse.elk.core.validation.IValidatingGraphElementVisitor;
import org.eclipse.elk.graph.properties.GraphFeature;
import org.eclipse.elk.graph.properties.IProperty;

/**
 * Data type used to store information for a layout algorithm.
 * 
 * @author msp
 */
public final class LayoutAlgorithmData implements ILayoutMetaData {

    /** identifier of the layout provider. */
    private final String id;
    /** user friendly name of the layout algorithm. */
    private final String name;
    /** runtime instance of the layout algorithm. */
    private final InstancePool<AbstractLayoutProvider> providerPool;
    /** layout category identifier. */
    private final String category;
    /** name of the bundle of which this algorithm is part of. */
    private final String melkBundleName;
    /** id of the (eclipse) bundle in which the melk file resides. */
    private final String definingBundleId;
    /** detail description. */
    private final String description;
    /** a path to a preview image for displaying in user interfaces. */
    private final String imagePath;
    /** Set of supported graph features. */
    private final Set<GraphFeature> supportedFeatures;
    /** Validator that can be applied to input graphs before the algorithm is executed. */
    private final Class<? extends IValidatingGraphElementVisitor> validatorClass;
    
    /** Map of known layout options. Keys are option IDs, values are the default values. */
    private final Map<String, Object> knownOptions = new HashMap<>();
    
    /**
     * Create a layout algorithm data entry.
     */
    private LayoutAlgorithmData(final Builder builder) {
        this.id = builder.id;
        this.name = builder.name;
        this.description = builder.description;
        this.providerPool = new InstancePool<AbstractLayoutProvider>(builder.providerFactory);
        this.category = builder.category;
        this.melkBundleName = builder.melkBundleName;
        this.definingBundleId = builder.definingBundleId;
        this.imagePath = builder.imagePath;
        if (builder.supportedFeatures == null) {
            this.supportedFeatures = EnumSet.noneOf(GraphFeature.class);
        } else {
            this.supportedFeatures = builder.supportedFeatures;
        }
        this.validatorClass = builder.validatorClass;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (obj instanceof LayoutAlgorithmData) {
            return this.id.equals(((LayoutAlgorithmData) obj).id);
        }
        return false;
    }
    
    @Override
    public int hashCode() {
        return id.hashCode();
    }

    @Override
    public String toString() {
        return "Layout Algorithm: " + id;
    }
    
    /**
     * Sets the knowledge status of the given layout option. This method should be called only from
     * {@link ILayoutMetaDataProvider} implementations.
     * 
     * @param property layout option
     * @param defaultValue the default value, or {@code null} if none is specified
     */
    public void addKnownOption(final IProperty<?> property, final Object defaultValue) {
        knownOptions.put(property.getId(), defaultValue);
    }
    
    /**
     * Returns the set of IDs of layout options declared to be known by this algorithm.
     * 
     * @return set of known layout option IDs.
     */
    public Set<String> getKnownOptionIds() {
        return knownOptions.keySet();
    }
    
    /**
     * Determines whether the layout algorithm knows the given layout option.
     * 
     * @param property layout option
     * @return true if the associated layout algorithm knows the option
     */
    public boolean knowsOption(final IProperty<?> property) {
        return knowsOption(property.getId());
    }
    
    /**
     * Determines whether the layout algorithm knows the given layout option.
     * 
     * @param propertyId layout option id.
     * @return true if the associated layout algorithm knows the option
     */
    public boolean knowsOption(final String propertyId) {
        return knownOptions.containsKey(propertyId);
    }
    
    /**
     * Returns the layout algorithm's default value for the given option.
     * 
     * @param property layout option
     * @return the associated default value, or {@code null} if there is none
     */
    public Object getDefaultValue(final IProperty<?> property) {
        return getDefaultValue(property.getId());
    }
    
    /**
     * Returns the layout algorithm's default value for the given option.
     * 
     * @param propertyId layout option id.
     * @return the associated default value, or {@code null} if there is none
     */
    public Object getDefaultValue(final String propertyId) {
        return knownOptions.get(propertyId);
    }
    
    /**
     * Check whether the given graph feature is supported.
     * 
     * @param graphFeature the graph feature to check
     * @return {@code true} if the algorithm supports the given feature
     */
    public boolean supportsFeature(final GraphFeature graphFeature) {
        return supportedFeatures.contains(graphFeature);
    }
    
    /**
     * @return the set of supported features
     */
    public Set<GraphFeature> getSupportedFeatures() {
        return supportedFeatures;
    }
    
    @Override
    public String getId() {
        return id;
    }

    @Override
    public String getName() {
        return name;
    }
    
    @Override
    public String getDescription() {
        return description;
    }

    /**
     * Return an instance pool for layout providers. If multiple threads execute the layout
     * algorithm in parallel, each thread should use its own instance of the algorithm.
     *
     * @return a layout provider instance pool
     */
    public InstancePool<AbstractLayoutProvider> getInstancePool() {
        return providerPool;
    }
    
    /**
     * Return a validator class that can be applied to input graphs before the algorithm is executed.
     * 
     * @return the graph validator class, or {@code null}
     */
    public Class<? extends IValidatingGraphElementVisitor> getValidatorClass() {
        return validatorClass;
    }

    /**
     * Returns the layout category identifier. Layout categories are represented by {@link LayoutCategoryData}.
     *
     * @return the category identifier
     */
    public String getCategoryId() {
        return category;
    }

    /**
     * Returns the bundle name. Bundles are used to group layout algorithms according to the library
     * in which they are contained.
     *
     * @return the bundle name
     */
    public String getBundleName() {
        return melkBundleName;
    }

    /**
     * @return the id of the (eclipse) bundle in which the melk file resides
     */
    public String getDefiningBundleId() {
        return definingBundleId;
    }
    
    /**
     * Returns the preview image path.
     * 
     * @return the preview image path
     */
    public String getPreviewImagePath() {
        return imagePath;
    }
    
    /**
     * Builder for {@link LayoutAlgorithmData} instances.
     */
    public static class Builder {

        private String id;
        private String name;
        private IFactory<AbstractLayoutProvider> providerFactory;
        private String category;
        private String melkBundleName;
        private String definingBundleId;
        private String description;
        private String imagePath;
        private Set<GraphFeature> supportedFeatures;
        private Class<? extends IValidatingGraphElementVisitor> validatorClass;
        
        /**
         * Create an instance with the configured values.
         */
        public LayoutAlgorithmData create() {
            return new LayoutAlgorithmData(this);
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getId() id}.
         */
        public Builder id(final String aid) {
            this.id = aid;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getName() name}.
         */
        public Builder name(final String aname) {
            this.name = aname;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getInstancePool() providerFactory}.
         */
        public Builder providerFactory(final IFactory<AbstractLayoutProvider> aproviderFactory) {
            this.providerFactory = aproviderFactory;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getCategoryId() category}.
         */
        public Builder category(final String acategory) {
            this.category = acategory;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getBundleName() melkBundleName}.
         */
        public Builder melkBundleName(final String amelkBundleName) {
            this.melkBundleName = amelkBundleName;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getDefiningBundleId() definingBundleId}.
         */
        public Builder definingBundleId(final String adefiningBundleId) {
            this.definingBundleId = adefiningBundleId;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getDescription() description}.
         */
        public Builder description(final String adescription) {
            this.description = adescription;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getPreviewImagePath() imagePath}.
         */
        public Builder imagePath(final String aimagePath) {
            this.imagePath = aimagePath;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getSupportedFeatures() supportedFeatures}.
         */
        public Builder supportedFeatures(final Set<GraphFeature> asupportedFeatures) {
            this.supportedFeatures = asupportedFeatures;
            return this;
        }
        
        /**
         * Configure the {@link LayoutAlgorithmData#getValidatorClass() validatorClass}.
         */
        public Builder validatorClass(final Class<? extends IValidatingGraphElementVisitor> avalidator) {
            this.validatorClass = avalidator;
            return this;
        }
        
    }

}
