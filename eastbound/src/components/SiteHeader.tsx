"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { useEffect, useState } from "react";

const links = [
  { href: "/plan", label: "Plan" },
  { href: "/destinations", label: "Destinations" },
  { href: "/itineraries", label: "Itineraries" },
  { href: "/transport", label: "Transport" },
  { href: "/essentials", label: "Essentials" },
];

export function SiteHeader() {
  const pathname = usePathname();
  const [scrolled, setScrolled] = useState(false);
  const [open, setOpen] = useState(false);
  const isHome = pathname === "/";

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 24);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  useEffect(() => {
    setOpen(false);
  }, [pathname]);

  const solid = !isHome || scrolled || open;

  return (
    <header className={`site-header${solid ? " is-solid" : ""}`}>
      <div className="site-header__inner">
        <Link href="/" className="brand" aria-label="EASTBOUND home">
          <span className="brand__mark">東</span>
          <span className="brand__name">EASTBOUND</span>
        </Link>

        <nav className="site-nav" aria-label="Primary">
          {links.map((link) => (
            <Link
              key={link.href}
              href={link.href}
              className={pathname.startsWith(link.href) ? "is-active" : undefined}
            >
              {link.label}
            </Link>
          ))}
        </nav>

        <Link href="/plan" className="header-cta">
          Start planning
        </Link>

        <button
          type="button"
          className="nav-toggle"
          aria-expanded={open}
          aria-controls="mobile-nav"
          onClick={() => setOpen((v) => !v)}
        >
          <span className="sr-only">Menu</span>
          <span />
          <span />
        </button>
      </div>

      <div
        id="mobile-nav"
        className={`mobile-nav${open ? " is-open" : ""}`}
        hidden={!open}
      >
        {links.map((link) => (
          <Link key={link.href} href={link.href}>
            {link.label}
          </Link>
        ))}
        <Link href="/plan" className="btn btn--primary">
          Start planning
        </Link>
      </div>
    </header>
  );
}
