/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.compaction.recthull;

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.ArrayList;


/**
 * Simple scaffold for a scanline algorithm.
 * 
 * @param <T>
 *            The type of the points the scanline works on.
 */
public final class Scanline<T> {

    /** The points to process along. */
    private List<T> points;
    /** The comparator to sort the points. */
    private Comparator<T> comparator;
    /** Handlers processing a point. */
    private List<EventHandler<T>> eventHandlers;

    private Scanline(final List<T> points, final Comparator<T> comparator,
            final Iterable<EventHandler<T>> eventHandlers) {
        this.comparator = comparator;
        this.points = points;
        this.eventHandlers = new ArrayList<>();
        eventHandlers.forEach(this.eventHandlers::add);
    }

    /**
     * @param points
     *            the points to process. They will be sorted according to {@link #comparator} and
     *            are then passed in the determined order to the passed {@link #eventHandlers}.
     * @param comparator
     *            a comparator to sort the passed points.
     * @param eventHandlers
     *            handlers to treat the points.
     * @param <T>
     *            The type of the points.
     */
    public static <T> void execute(final Iterable<T> points, final Comparator<T> comparator,
            final Iterable<EventHandler<T>> eventHandlers) {
        // copy the points! we will resort them!
        List<T> copy = new ArrayList<>();
        points.forEach(copy::add);
        new Scanline<T>(copy, comparator, eventHandlers).go();
    }

    /**
     * @param points
     *            the points to process. They will be sorted according to {@link #comparator} and
     *            are then passed in the determined order to the passed {@link #eventHandlers}.
     * @param comparator
     *            a comparator to sort the passed points.
     * @param eventHandler
     *            a handler to treat the points.
     * @param <T>
     *            The type of the point.
     */
    public static <T> void execute(final Iterable<T> points, final Comparator<T> comparator,
            final EventHandler<T> eventHandler) {
        execute(points, comparator, Arrays.asList(eventHandler));
    }

    private void go() {

        // sort
        Collections.sort(points, comparator);

        // now move scanline
        for (T p : points) {
            for (EventHandler<T> h : eventHandlers) {
                h.handle(p);
            }
        }

    }

    /**
     * An event handler, gets passed a point of type T and does with it whatever it likes.
     * 
     * @author uru
     * @param <T>
     *            The type of the point.
     */
    @FunctionalInterface
    public interface EventHandler<T> {
        /**
         * @param p
         *            the latest point the scanline algorithm found.
         */
        void handle(T p);
    }

}
