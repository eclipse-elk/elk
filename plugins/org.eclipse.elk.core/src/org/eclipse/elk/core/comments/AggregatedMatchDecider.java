/*******************************************************************************
 * Copyright (c) 2016, 2017 Kiel University and others.
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
package org.eclipse.elk.core.comments;

import java.util.Collection;
import java.util.Collections;
import java.util.Map;
import java.util.Objects;
import java.util.function.ToDoubleFunction;
import java.util.stream.Collectors;

/**
 * Selects the attachment target with the best aggregated heuristic result. The matcher results are first aggregated
 * using a configurable aggregation function. Then, the attachment target with the highest result is selected, provided
 * that it is at least as high as or higher than a configurable lower bound. The class provides a number of pre-defined
 * aggregation functions.
 * 
 * @param <T>
 *            type of attachment targets.
 */
public final class AggregatedMatchDecider<T> implements IDecider<T> {
    
    /** Aggregator used to aggregate heuristics results. */
    private ToDoubleFunction<Collection<Double>> aggregator = AggregatedMatchDecider::max;
    /** The minimum aggregate result for a comment to be attached to anything. */
    private double lowerBoundary = 0.0;
    /** Whether an attachment is accepted if the aggregate value is exactly the lower boundary. */
    private boolean includeLowerBoundary = false;
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /**
     * Configures the attachment decider to use the given function to aggregate heuristic results.
     * 
     * <p>
     * If this method is not called, the {@link #max(Collection) maximum aggregator} will be used.
     * </p>
     * 
     * @param f
     *            the aggregator to use. There are no constraints on the values it produces, except that they must not
     *            be negative.
     * @return this object for method chaining.
     * @throws IllegalArgumentException
     *             if the function is {@code null}.
     */
    public AggregatedMatchDecider<T> withAggregator(final ToDoubleFunction<Collection<Double>> f) {
        Objects.requireNonNull(f, "Aggregator cannot be null.");
        
        aggregator = f;
        return this;
    }
    
    /**
     * Configures the attachment decider to impose the given lower boundary on attachability decisions. This means that
     * if the attachment target with the highest aggregated result doesn't score at least as high as the lower boundary,
     * it is not returned.
     * 
     * <p>
     * If this method is not called, all attachment targets will be considered eligible.
     * </p>
     * 
     * @param lower
     *            the lower boundary, {@code >= 0}.
     * @return this object for method chaining.
     * @throws IllegalArgumentException
     *             if the lower boundary is {@code < 0}.
     */
    public AggregatedMatchDecider<T> withLowerAttachmentBoundary(final double lower) {
        if (lower < 0) {
            throw new IllegalArgumentException("Lower boundary must be >= 0.");
        }
        
        lowerBoundary = lower;
        
        return this;
    }
    
    /**
     * Configures whether the attachment decider allows attachments already if the highest aggregated result is exactly
     * the lower boundary.
     * 
     * <p>
     * If this method is not called, the decider requires the highest aggregated result to be strictly greater than the
     * lower boundary.
     * </p>
     * 
     * @param include
     *            {@code true} if the lower boundary should be included.
     * @return this object for method chaining.
     */
    public AggregatedMatchDecider<T> withLowerBoundaryIncluded(final boolean include) {
        includeLowerBoundary = include;
        return this;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // IAttachmentDecider
    
    @Override
    public T makeAttachmentDecision(final Map<T, Map<Class<? extends IMatcher<?, T>>, Double>> normalizedHeuristics) {
        double max = Double.NEGATIVE_INFINITY;
        T maxTarget = null;
        
        for (Map.Entry<T, Map<Class<? extends IMatcher<?, T>>, Double>> entry : normalizedHeuristics.entrySet()) {
            double aggregate = aggregator.applyAsDouble(entry.getValue().values());
            if (aggregate < 0) {
                throw new IllegalStateException("The aggregator provided a value < 0.");
            }
            
            if (aggregate > max) {
                max = aggregate;
                maxTarget = entry.getKey();
            }
        }
        
        if (includeLowerBoundary) {
            return max >= lowerBoundary ? maxTarget : null;
        } else {
            return max > lowerBoundary ? maxTarget : null;
        }
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////
    // Aggregation Functions
    
    /**
     * Aggregator that returns the maximum heuristic result.
     * 
     * @param values heuristic results.
     * @return the maximum result.
     */
    public static double max(final Collection<Double> values) {
        return Collections.max(values);
    }
    
    /**
     * Aggregator that returns the minimum heuristic result.
     * 
     * @param values heuristic results.
     * @return the minimum result.
     */
    public static double min(final Collection<Double> values) {
        return Collections.min(values);
    }
    
    /**
     * Aggregator that returns the average over all heuristic results.
     * 
     * @param values heuristic results.
     * @return the average.
     */
    public static double avg(final Collection<Double> values) {
        return values.stream().collect(Collectors.averagingDouble((d) -> d));
    }
    
    /**
     * Aggregator that returns the sum of all heuristic results.
     * 
     * @param values heuristic results.
     * @return the sum.
     */
    public static double sum(final Collection<Double> values) {
        return values.stream().collect(Collectors.summingDouble((d) -> d));
    }

}
