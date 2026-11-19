/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
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
package org.eclipse.elk.core.options;

import java.util.EnumSet;

/**
 * Determines how the content of compound nodes is to be aligned if the compound node's size exceeds
 * the bounding box of the content (i.e. child nodes). This might be the case if for a compound node
 * the {@link org.eclipse.elk.core.options.SizeConstraint SizeConstraint} of
 * {@link org.eclipse.elk.core.options.SizeConstraint#MINIMUM_SIZE MINIMUM_SIZE} is set and
 * {@link org.eclipse.elk.core.options.LayeredOptions#MIN_WIDTH MIN_WIDTH} and
 * {@link org.eclipse.elk.core.options.LayeredOptions#MIN_HEIGHT MIN_HEIGHT} are set large enough.
 * 
 * <p>This property is to be used as an {@link java.util.EnumSet EnumSet}; it should be comprised of
 * exactly one value prefixed with {@code V_} and one prefixed with {@code H_}.</p>
 * 
 * <p>Default values are {@link ContentAlignment#V_TOP V_TOP} and {@link ContentAlignment#H_LEFT
 * H_LEFT}.</p>
 * 
 * @author uru
 */
public enum ContentAlignment {

    /** Content should be vertically aligned to the top. */
    V_TOP,
    /** Content should be vertically centered. */
    V_CENTER,
    /** Content should be vertically aligned to the bottom. */
    V_BOTTOM,
    /** Content should be horizontally aligned to the left. */
    H_LEFT,
    /** Content should be horizontally centered. */
    H_CENTER,
    /** Content should be horizontally aligned to the right. */
    H_RIGHT;

    /**
     * @return a set containing {@link #V_CENTER} and {@link #H_CENTER}.
     */
    public static EnumSet<ContentAlignment> centerCenter() {
        return EnumSet.of(V_CENTER, H_CENTER);
    }

    /**
     * @return a set containing {@link #V_TOP} and {@link #H_LEFT}.
     */
    public static EnumSet<ContentAlignment> topLeft() {
        return EnumSet.of(V_TOP, H_LEFT);
    }

    /**
     * @return a set containing {@link #V_BOTTOM} and {@link #H_RIGHT}.
     */
    public static EnumSet<ContentAlignment> bottomRight() {
        return EnumSet.of(V_BOTTOM, H_RIGHT);
    }

    /**
     * @return a set containing {@link #V_TOP} and {@link #H_CENTER}.
     */
    public static EnumSet<ContentAlignment> topCenter() {
        return EnumSet.of(V_TOP, H_CENTER);
    }
}
