/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.comments;

import java.awt.geom.Rectangle2D;
import java.util.Objects;

/**
 * Determines if a comment is eligible for attachment based on its size. Use the methods named {@code withXXX} to
 * configure the filter.
 * 
 * @param <C>
 *            type of comments.
 */
public class SizeFilter<C> implements IFilter<C> {
    
    /** The bounds provider to use. */
    private IBoundsProvider<C, ?> boundsProvider = null;
    /** The maximum area for a comment to still be eligible. */
    private double maxArea = -1;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /**
     * Configures the filter to consider comments up to the given are to be attachable.
     * 
     * <p>
     * If this method is not called, no comment is considered eligible for attachment.
     * </p>
     * 
     * @param area
     *            the maximum possible area.
     * @return this object for method chaining.
     */
    public SizeFilter<C> withMaximumArea(final double area) {
        if (area < 0) {
            throw new IllegalArgumentException("Maximum area must be >= 0.");
        }
        
        maxArea = area;
        return this;
    }
    
    /**
     * Configures the filter to use the given bounds provider to determine the bounds of comments.
     * 
     * <p>
     * If this method is not called, the filter will throw an exception during preprocessing.
     * </p>
     * 
     * @param provider
     *            the bounds provider to use.
     * @return this object for method chaining.
     */
    public SizeFilter<C> withBoundsProvider(final IBoundsProvider<C, ?> provider) {
        Objects.requireNonNull(provider, "Bounds provider must not be null.");
        
        this.boundsProvider = provider;
        return this;
    }
    
    /**
     * Checks whether the current configuration is valid.
     * 
     * @throws IllegalStateException
     *             if the configuration is invalid.
     */
    private void checkConfiguration() {
        if (boundsProvider == null) {
            throw new IllegalStateException("A bounds provider is required.");
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // IFilter
    
    @Override
    public void preprocess(final IDataProvider<C, ?> dataProvider, final boolean includeHierarchy) {
        checkConfiguration();
    }

    @Override
    public boolean eligibleForAttachment(final C comment) {
        Rectangle2D.Double bounds = boundsProvider.boundsForComment(comment);
        return bounds.height * bounds.width <= maxArea;
    }

}
