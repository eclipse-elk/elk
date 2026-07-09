/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.mrtree.options;

/**
 * 
 * @author jnc, sdo
 */
public enum CompactionMode {
    /** No Compaction. */
    NONE,
    /** One dimensional Compaction that preserves the spacings between levels. */
    LEVEL_PRESERVING,
    /** One dimensional Compaction that does not preserve the spacings between levels. 
     * WARNING: Do not use with EdgeRoutingMode AvoidOverlap */
    AGGRESSIVE
}
